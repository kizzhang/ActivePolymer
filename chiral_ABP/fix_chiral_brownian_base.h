/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_CHIRAL_BROWNIAN_BASE_H
#define LMP_FIX_CHIRAL_BROWNIAN_BASE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChiralBrownianBase : public Fix {
 public:
  FixChiralBrownianBase(class LAMMPS *, int, char **);
  ~FixChiralBrownianBase() override;
  void init() override;
  int setmask() override;
  void reset_dt() override;

 protected:
  int seed;                  // RNG seed
  double dt, sqrtdt;         // time step interval and its sqrt
  int gamma_t_flag;          // 0/1 if isotropic translational damping is unset/set
  int gamma_r_flag;          // 0/1 if isotropic rotational damping is unset/set
  int rot_temp_flag;         // 0/1 if rotational temperature is unset/set
  int planar_rot_flag;       // 0/1 if rotation is constrained to 2D (xy) plane

  double gamma_t, gamma_r;    // translational and rotational (isotropic) damping params

  int gaussian_noise_flag;    // 0/1 for uniform/gaussian noise
  double *gaussian_mean;      // mean of gaussian noise

  double temp;        // temperature
  double rot_temp;    // temperature
  double g1, g2;      // prefactors in time stepping

  class RanMars *rng;
};

}    // namespace LAMMPS_NS
#endif
