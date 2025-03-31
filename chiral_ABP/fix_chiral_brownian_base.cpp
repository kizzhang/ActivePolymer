/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified from fix_brownian_base.cpp.

   Contributing author: Zhiyu Zhang (City University of Hong Kong)
------------------------------------------------------------------------- */

#include "fix_chiral_brownian_base.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
FixChiralBrownianBase::FixChiralBrownianBase(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), gaussian_mean(nullptr), rng(nullptr)
{
  time_integrate = 1;

  rot_temp_flag = 0;
  planar_rot_flag = 0;
  gamma_t_flag = gamma_r_flag = 0;
  rot_temp_flag = 0;
  g2 = 0.0;

  if (narg < 9) utils::missing_cmd_args(FLERR, "fix chiral brownian", error);

  temp = utils::numeric(FLERR, arg[3], false, lmp);
  if (temp <= 0) error->all(FLERR, "Fix chiral brownian temp must be > 0.0");

  seed = utils::inumeric(FLERR, arg[4], false, lmp);
  if (seed <= 0) error->all(FLERR, "Fix chiral brownian seed must be > 0");
  
  int iarg = 5;
  if (strcmp(arg[iarg], "chiral") == 0) {
    if (narg < iarg + 3) utils::missing_cmd_args(FLERR, "fix chiral brownian gaussian mean", error);

    delete[] gaussian_mean;
    gaussian_mean = new double[3];

    gaussian_mean[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    gaussian_mean[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
    gaussian_mean[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
    iarg = iarg + 4;
  } 
  else error->all(FLERR, "Expected 'chiral' keyword for fix chiral brownian");

  while (iarg < narg) {
    if (strcmp(arg[iarg], "gamma_t") == 0) {
      if (narg == iarg + 1) { error->all(FLERR, "Illegal fix chiral brownian command."); }

      gamma_t_flag = 1;
      gamma_t = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (gamma_t <= 0) error->all(FLERR, "Fix chiral brownian gamma_t must be > 0.");
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "gamma_r") == 0) {
      if (narg == iarg + 1) { error->all(FLERR, "Illegal fix chiral brownian command."); }

      gamma_r_flag = 1;
      gamma_r = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (gamma_r <= 0) error->all(FLERR, "Fix chiral brownian gamma_r must be > 0.");
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "rotation_temp") == 0) {
      if (narg == iarg + 1) { error->all(FLERR, "Illegal fix chiral brownian command."); }

      rot_temp_flag = 1;
      rot_temp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (rot_temp <= 0) error->all(FLERR, "Fix chiral brownian rotation_temp must be > 0.");
      iarg = iarg + 2;

    } else if (strcmp(arg[iarg], "planar_rotation") == 0) {

      planar_rot_flag = 1;
      if (domain->dimension == 2)
        error->all(FLERR, "The planar_rotation keyword is not allowed for 2D simulations");
      iarg = iarg + 1;

    } else {
      error->all(FLERR, "Illegal fix chiral brownian command.");
    }
  }
  if (!rot_temp_flag) rot_temp = temp;

  // initialize Marsaglia RNG with processor-unique seed
  rng = new RanMars(lmp, seed + comm->me);
}

/* ---------------------------------------------------------------------- */

int FixChiralBrownianBase::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixChiralBrownianBase::~FixChiralBrownianBase()
{
  delete[] gaussian_mean;
  delete rng;
}

/* ---------------------------------------------------------------------- */

void FixChiralBrownianBase::init()
{
  dt = update->dt;
  sqrtdt = sqrt(dt);
  g1 = force->ftm2v;
  g2 = sqrt(2 * force->boltz / dt / force->mvv2e);
}

void FixChiralBrownianBase::reset_dt()
{
  double sqrtdt_old = sqrtdt;
  dt = update->dt;
  sqrtdt = sqrt(dt);
  g2 *= sqrtdt_old / sqrtdt;
}
