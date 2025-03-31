#ifdef FIX_CLASS

FixStyle(tangential/propel, FixTangentialPropel)

#else

#ifndef LMP_FIX_TANGENTIAL_PROPEL_H
#define LMP_FIX_TANGENTIAL_PROPEL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTangentialPropel : public Fix {
 public:
  FixTangentialPropel(class LAMMPS *, int, char **);

  void init() override;
  int setmask() override;
  void setup(int) override;
  void post_force(int) override;

 private:
  double magnitude;   // Active force magnitude
  int angletype;      // Angle type to consider for tangential force
  int mode;           // Mode of operation (TANGENT)
  
  void post_force_tangent(int);
};

} // namespace LAMMPS_NS

#endif
#endif