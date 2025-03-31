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

#include "fix_tangential_propel.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;


static constexpr double TOL = 1e-14;

/* ---------------------------------------------------------------------- */

FixTangentialPropel::FixTangentialPropel(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  virial_global_flag = virial_peratom_flag = 1;

  if (narg != 5) error->all(FLERR, "Illegal fix tangential/propel command: Syntax: fix ID group tangential/propel magnitude angle_type");

  //thermo_virial = 1;
  
  magnitude = utils::numeric(FLERR, arg[3], false, lmp);
  angletype = utils::inumeric(FLERR, arg[4], false, lmp);
}

/* ---------------------------------------------------------------------- */

int FixTangentialPropel::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTangentialPropel::init() {
  if (!neighbor->anglelist)
      error->all(FLERR, "Fix tangential/propel tangent requires angle interactions"); 
  if (angletype < 1 || angletype > atom->nangletypes)
      error->all(FLERR, "Invalid angle type for fix tangential/propel tangent");
  
}

/* ---------------------------------------------------------------------- */

void FixTangentialPropel::setup(int vflag) {
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTangentialPropel::post_force(int vflag)
{
  post_force_tangent(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTangentialPropel::post_force_tangent(int vflag) {
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  imageint *image = atom->image;
  double unwrap[3];
  double vi[6];

  if (vflag) v_setup(vflag);
  else evflag = 0;

  for (int n = 0; n < nanglelist; ++n) {
    int i = anglelist[n][0];
    int j = anglelist[n][1];
    int k = anglelist[n][2];
    int type = anglelist[n][3];

    if (type != angletype) continue;
    if (!(mask[j] & groupbit)) continue;

    // Compute tangent vector from i to k
    double delx = x[k][0] - x[i][0];
    double dely = x[k][1] - x[i][1];
    double delz = x[k][2] - x[i][2];

    double rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < TOL) continue;

    double rinv = 1.0 / sqrt(rsq);
    delx *= rinv;
    dely *= rinv;
    delz *= rinv;

    // Apply force to atom j
    double fx = magnitude * delx;
    double fy = magnitude * dely;
    double fz = magnitude * delz;

    f[j][0] += fx;
    f[j][1] += fy;
    f[j][2] += fz;

    // Virial contribution
    if (evflag) {
      domain->unmap(x[j], image[j], unwrap);
      vi[0] = fx * unwrap[0];
      vi[1] = fy * unwrap[1];
      vi[2] = fz * unwrap[2];
      vi[3] = fx * unwrap[1];
      vi[4] = fx * unwrap[2];
      vi[5] = fy * unwrap[2];
      v_tally(j, vi);
    }
  }
}