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
   Originally modified from fix_brownian_sphere.cpp.

   Contributing author: Zhiyu Zhang (City University of Hong Kong)
------------------------------------------------------------------------- */

#include "fix_chiral_brownian_sphere.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChiralBrownianSphere::FixChiralBrownianSphere(LAMMPS *lmp, int narg, char **arg) :
    FixChiralBrownianBase(lmp, narg, arg)
{
  if (!gamma_t_flag || !gamma_r_flag) error->all(FLERR, "Illegal fix chiral/brownian/sphere command.");
  if (!atom->mu_flag) error->all(FLERR, "Fix chiral/brownian/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixChiralBrownianSphere::init()
{
  FixChiralBrownianBase::init();

  g3 = g1 / gamma_r;
  g4 = g2 * sqrt(rot_temp / gamma_r);
  g1 /= gamma_t;
  g2 *= sqrt(temp / gamma_t);
}

/* ---------------------------------------------------------------------- */

void FixChiralBrownianSphere::initial_integrate(int /*vflag */)
{
  if (domain->dimension == 2) {

    initial_integrate_templated<1, 0>();
    
  } else if (planar_rot_flag) {

    initial_integrate_templated<0, 1>();
    
  } else {

    initial_integrate_templated<0, 0>();
    
  }
}

/* ---------------------------------------------------------------------- */

template <int Tp_2D, int Tp_2Drot>
void FixChiralBrownianSphere::initial_integrate_templated()
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double wx, wy, wz;
  double **torque = atom->torque;
  double **mu = atom->mu;
  double mux, muy, muz, mulen;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dx, dy, dz;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (Tp_2D) {            // 2D ROTATION
        dz = 0;
        wx = wy = 0;
        
        dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
        dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
        wz = rng->gaussian(gaussian_mean[2],1) * g4;

      } else if (Tp_2Drot) {  // 3D PLANAR ROTATION
        wx = wy = 0;
        
        dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
        dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
        dz = dt * (g1 * f[i][2] + g2 * rng->gaussian());
        wz = rng->gaussian(gaussian_mean[2],1) * g4;

      } else {                // 3D ROTATION
        dx = dt * (g1 * f[i][0] + g2 * rng->gaussian());
        dy = dt * (g1 * f[i][1] + g2 * rng->gaussian());
        dz = dt * (g1 * f[i][2] + g2 * rng->gaussian());
        wx = rng->gaussian(gaussian_mean[0],1) * g4;
        wy = rng->gaussian(gaussian_mean[1],1) * g4;
        wz = rng->gaussian(gaussian_mean[2],1) * g4;
      }
    
      x[i][0] += dx;
      v[i][0] = dx / dt;

      x[i][1] += dy;
      v[i][1] = dy / dt;

      x[i][2] += dz;
      v[i][2] = dz / dt;

      wx += g3 * torque[i][0];
      wy += g3 * torque[i][1];
      wz += g3 * torque[i][2];

      // store length of dipole as we need to convert it to a unit vector and
      // then back again

      mulen = sqrt(mu[i][0] * mu[i][0] + mu[i][1] * mu[i][1] + mu[i][2] * mu[i][2]);

      // unit vector at time t
      mux = mu[i][0] / mulen;
      muy = mu[i][1] / mulen;
      muz = mu[i][2] / mulen;

      // un-normalised unit vector at time t + dt
      mu[i][0] = mux + (wy * muz - wz * muy) * dt;
      mu[i][1] = muy + (wz * mux - wx * muz) * dt;
      mu[i][2] = muz + (wx * muy - wy * mux) * dt;

      // normalisation introduces the stochastic drift term due to changing from
      // Stratonovich to Ito interpretation
      MathExtra::norm3(mu[i]);

      // multiply by original magnitude to obtain dipole of same length
      mu[i][0] = mu[i][0] * mulen;
      mu[i][1] = mu[i][1] * mulen;
      mu[i][2] = mu[i][2] * mulen;
    }
  }
}
