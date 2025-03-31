# Active Polymer Simulation

## 1. LAMMPS source code
### 1.1. Tangentially-driven active polymer
Tangentially-driven by active force, active polymers are prone to collapse in closed topology (unknotted ring). ![image](https://github.com/user-attachments/assets/f9ae74b8-3651-4d3f-8814-2ec23f109fbb "Active force is along the tangent vector.")

The code is modified from the fix_propel_self.cpp in [LAMMPS](https://github.com/lammps).

### 1.2. Achiral active Brownian polymer 
Active Brownian particles can be chained together as a polymer. We can then run the simulation using fix/self/propel in LAMMPS. I will attach an **in.run** script for this later. Check my HOOMD-Blue code later.

### 1.3. Chiral active Brownian polymer
Chirality can be introduced into active Brownian particles. See this inclusive [paper](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.062120) by Sevilla. Chiral ABPs can also be chained to form a polymer. I will later upload LAMMPS source codes after benchmarking.


## 2. HOOMD scripts


