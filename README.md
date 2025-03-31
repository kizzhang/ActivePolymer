# Active Polymer Simulation
<mark>**EMBRACE OPEN-ACCESS SCIENCE!**</mark>

##  LAMMPS Simulation
#### • Tangentially-driven active polymer
Tangentially-driven by active force, active polymers are prone to [collapse](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.097801) (E. Locatelli et al.2021 PRL) in closed topology (unknotted ring).  

![image](https://github.com/user-attachments/assets/f9ae74b8-3651-4d3f-8814-2ec23f109fbb "Active force is along the tangent vector.")

The code is modified from the fix_propel_self.cpp in [LAMMPS](https://github.com/lammps).

#### • Achiral active Brownian polymer 
Active Brownian particles can be chained together as a polymer. See image below.

![active Brownian Particle](https://github.com/user-attachments/assets/b4a1933e-60c6-443c-add8-8b685f8075b2)

We can then run the simulation using `fix self/propel` in LAMMPS. I will post HOOMD-Blue code later.

#### • Chiral active Brownian polymer
Chirality can be introduced into active Brownian particles. See this comprehensive [paper](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.062120) by Sevilla. Chiral ABPs can also be chained to form a polymer. 

![chiral](https://github.com/user-attachments/assets/1825c94c-60d8-4b5c-8226-46b5ca1158da "Chiral active polymer with torque acting to rotate the direction of active force.")

As stated in the paper by Sevilla, the chirality is introduced by making the mean of rotational Gaussian noise some non-zero constant.

## Usage 
To compile the source codes, simply put the **.cpp file and its .h header file** into your `lammps/src` directory, and compile as detailed [here](https://docs.lammps.org/Build.html).

---
#### • Tangentially-driven active polymer 

Use 

```
fix         tang all tangential/propel ${mag} ${angle_tyle}
fix         brownian all brownian/sphere ${temp} ${seed} gamma_t ${gamma_t} gamma_r ${gamma_r}
```

**OR** for Langevin simulation:

```
fix         tang all tangential/propel ${mag} ${angle_tyle}
fix         langevin all langevin ${temp} ${temp} ${damp} ${seed}
```


(Active force is defined with `atom->anglelist` in LAMMPS). 

To run tutorial code: go to `tangential/run_scripts`.

---
#### • Achiral active Brownian polymer

Use 
```
fix         brownian all brownian/sphere ${temp} 4928459 gamma_t ${gamma_t} gamma_r ${gamma_r}
fix         actforce all propel/self dipole ${f_mag}
```

To run tutorial code: go to 'achiral/'

---
#### • Chiral active Brownian polymer

Use 

```
fix         chiral all chiral/brownian/sphere ${temp} ${seed} chiral ${tau_x} ${tau_y} ${tau_z} gamma_t ${gamma_t} gamma_r ${gamma_r}
fix         actforce all propel/self dipole ${f_mag}
```

`${tau_x,y,z}` is the torque magnitude along x,y,z axis. It also is the *mean* value for rotational Gaussian noise. You can use `plannar_rotation` flag to introduce chirality only in x-y plane.

To run tutorial code: go to 'chiral/run_scripts'

----
### NOTE: **`fix tangential/propel` does not require dipole. Fix with `dipole` should be used in parallel with `fix brownian/sphere`. Check LAMMPS for more details.**


## HOOMD scripts


## Citation
You can cite this repository as 
`Zhiyu Z., Active polymer simulation, GitHub repository 2025, https://github.com/kizzhang/ActivePolymer`
```
@misc{Zhiyu2025,
  author = {Zhiyu, Z.},
  title = {Active Polymer Simulation},
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/kizzhang/ActivePolymer}},
  commit = {main}
}
```

