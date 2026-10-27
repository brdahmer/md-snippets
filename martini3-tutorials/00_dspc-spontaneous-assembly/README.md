Reproduction of [Lipid bilayers - Part I – Self-assembly](https://cgmartini.nl/docs/tutorials/Martini3/LipidsI/) tutorial for Martini

This tutorial covers:
1. Spontaneous self-assembly of a distearoyl-phosphatidylcholine (DSPC) bilayer and a quick assessment of various properties that characterize lipid bilayers
2. How changing the nature of lipid head groups and tails influences these properties

As per the original tutorial, all Martini3 Lipid Parameter files can be downloaded from [here](https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters). This tutorial only uses `martini_v3.0.0_ffbonded_v2.itp` and `martini_v3.0.0_phospholipids_PC_v2.itp` besides the [Martini3 release](https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip). I like to keep all definitions within the same folder, so move the individual `.itp` files to the `martini_v300/` directory
___
### Bilayer self-assembly

1. Create a simulation box containing a random configuration of 128 DSPC lipids from the `.gro` file of a single (energy minimized) DSPC molecule:

```
gmx insert-molecules -ci DSPC-em.gro -box 7.5 7.5 7.5 -nmol 128 -radius 0.21 -try 500 -o 128_noW.gro
```

With this command, the `insert-molecules` module attempts to randomly place a given molecule (`-ci DSPC-em.gro`) multiple times (`-nmol` $*$ `-try`) within a simulation box (`-box 7.5 7.5 7.5`, \[x y z]).  The van der Waals radii (`-radius`) has to be increased from its default atomistic length (0.105 nm) to a value reflecting the size of Martini CG beads (0.21 nm).

2. Add 6 CG water molecules per lipid (768) from a `.gro` file of an equilibrated water box serving as template:

```
gmx solvate -cp 128_noW.gro -cs water.gro -o waterbox.gro -maxsol 768 -radius 0.21
```

3. Simple energy minimization of the solvated system:

First, build the `topology` file (`.top`):
```
#include "./martini_v300/martini_v3.0.0.itp" ; the particle definitions should be included first
#include "./martini_v300/martini_v3.0.0_ffbonded_v2.itp" ; the general lipid bonded definitions should be included before the specific lipid topologies
#include "./martini_v300/martini_v3.0.0_phospholipids_PC_v2.itp" ; remaining relevant .itps for molecules
#include "./martini_v300/martini_v3.0.0_solvents_v1.itp" ; for water

[ system ]
; Title
DSPC BILAYER SELF-ASSEMBLY

[ molecules ]
; Molecule types and their numbers are written in the order
; they appear in the structure file
DSPC 128 ; from gmx insert-molecules
W    768 ; from gmx solvate
```

And the `.mdp` file:
```
integrator               = steep
nsteps                   = 10000

nstxout                  = 0
nstvout                  = 0
nstfout                  = 0

cutoff-scheme            = Verlet
nstlist                  = 20
nsttcouple               = 20
nstpcouple               = 20
rlist                    = 1.35
verlet-buffer-tolerance  = -1
ns_type                  = grid
pbc                      = xyz

coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15	; 2.5 (with polarizable water)
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

constraints              = none
constraint_algorithm     = Lincs
lincs_order              = 8
lincs_warnangle          = 90
lincs_iter               = 2
```
Set up and run:
```
gmx grompp -f minimization.mdp -c waterbox.gro -p dspc.top -o em.tpr
gmx mdrun -s em.tpr -v -c em.gro
```

4. Self-assembly MD simulation:

```
gmx grompp -f martini_md.mdp -c em.gro -p dspc.top -o dspc-md.tpr
gmx mdrun -s dspc-md.tpr -v -x dspc-md.xtc -c dspc-md.gro
```

5. For visualization using PyMOL, it is recommended that the trajectory be converted to `.pdb` format because Martini bonds are much longer than usual, which messes up visualization. Select `System` OR `DSPC`

```
echo "0" | gmx trjconv -s dspc-md.tpr -f dspc-md.xtc -o dspc-md.pdb -pbc whole -conect
pymol dspc-md.pdb 
```

#### Important
If a bilayer did form during the simulation, check if it is normal to the **z-axis** (i.e., that the membrane is situated in the **xy-plane**). 
This can be done by profiling the `density` along all three axes; the **normal** is found in whatever plot shows high lipid content where water is low (water-free region in a hydrophobic core)

`echo "2" "3" | gmx density -s dspc-md.tpr -f dspc-md.gro -d X -ng 2 -o dens_x.xvg`
`echo "2" "3" | gmx density -s dspc-md.tpr -f dspc-md.gro -d Y -ng 2 -o dens_y.xvg`
`echo "2" "3" | gmx density -s dspc-md.tpr -f dspc-md.gro -d Z -ng 2 -o dens_z.xvg`

In our case, the bilayer was formed normal to the **x-axis!**

In cases like this, we can rotate the system with `gmx editconf -rotate`. 
For a normal along **x**, we rotate 90 degrees along **y**. 
`gmx editconf -f dspc-md.gro -o dspc-rot.gro -rotate 0 90 0`
For a normal along **y**, we rotate 90 degrees along **x**:
`gmx editconf -f dspc-md.gro -o dspc-rot.gro -rotate 90 0 0`
For the trajectory:
`gmx trjconv -s dspc-md.tpr -f dspc-rot.gro -o dspc-rot-wrapped.gro -pbc mol -ur compact`

>Verify the orientation:

`gmx density -s dspc-md.tpr -f dspc-rot-wrapped.gro -d Z -ng 2 -o dens_z_rot.xvg`


### Bilayer equilibration and analysis

```
mkdir equilibration
cp dspc-rot-wrapped.gro martini_md.mdp dspc.top equilibration/
```

1. The above command copies the necessary files to a new directory:
	1. The (corrected) output of the self-assembly simulation
	2. The `dspc.top` topology file (modify the `#include` statements at the top to point to the right files)
	3. The `martini_md.mdp` file. Change: 
		1. `comm-grps = DSPC W`
		2. `Pcoupltype = semi-isotropic` (We'll keep the C-rescale barostat with a 4 ps `tau_p`. Parrinello-Rahman is also suitable, with `tau_p = 12`)
		3. `tau_t = 1.0 1.0
		4. `ref_t = 340 340`
		5. `ref_p = 1.0 1.0`
		6. `compressibility = 3e-4 3e-4`
		7. `gen_vel = yes` 
		   `gen_temp = 340`
		   `gen_seed = -1`
		   to
		   `gen_vel = no`
		   `continuation = yes`

2. Set up and run:
```
gmx grompp -f martini_md.mdp -c dspc-rot-wrapped.gro -p dspc.top -o dspc-eq.tpr
gmx mdrun -s dspc-eq.tpr -v -x dspc-eq.xtc -c dspc-eq.gro
```
___
#### Area per lipid

Generally, the area per lipid (APL) can be calculated by dividing the area of the simulation box by half the number of lipids in the system. The box dimensions (along the x and y axes) can be obtained by running `gmx energy`, selecting the Box-X and Box-Y options:
`gmx energy -f ener.edr -o box-xy.xvg`
`11 12 0`
`qtgrace -nxy box-xy.xvg`

<img width="618" height="565" alt="Pasted image 20261027150947" src="https://github.com/user-attachments/assets/ec3c7062-cafc-4037-a7e5-12dc17082231" />

Because of the semi-isotropic coupling, both dimensions vary simultaneously and thus their output is the same. Looking at the graph, we can discard the first 15000 ps from the simulation as non-equilibrated (relaxation time).

To calculate APL, run the following in a `python3` session. This will create an `area.xvg` file.

```
import numpy as np
xy = np.loadtxt('box-xy.xvg', comments=('#','@'))
areas = xy[:,1] * xy[:,2] / 64
np.savetxt('area.xvg', np.column_stack((xy[:,0], areas)))
quit()
```

<img width="793" height="552" alt="Pasted image 20261027151322" src="https://github.com/user-attachments/assets/d9992761-9174-4c31-a37d-7d107c4301d4" />

The averages and error estimates can be extracted with the `gmx analyze` module:

`
`gmx analyze -f area.xvg -ee`
Or, discarding the first 15000 ps:
`gmx analyze -f area.xvg -ee -b 15000`

<img width="687" height="620" alt="Pasted image 20261027151958" src="https://github.com/user-attachments/assets/cd94f006-59c8-4803-9208-5392dca727db" />

This gives APL of 0.593 ± 0.003 nm². Experimental data using small-angle neutron scattering give an APL of 0.65 ± 0.05 nm². Given the short simulation time, small hydration, and possible assymetry in lipid numbers, this is fine. 
___
#### Bilayer thickness

Bilayer thickness can be extracted using `gmx density`. First, however, we'll create an index (`.ndx`) file for whichever functional groups in the lipid we want to use (phosphate and ammonium headgroup beads, carbon tail beads, etc.). For phosphate:

`gmx make_ndx -f dspc-eq.gro`
```
a P*
q
```
`echo "4" | gmx density -f dspc-eq.xtc -s dspc-eq.tpr -b 15000 -n index.ndx -o p-density.xvg -xvg no`

<img width="759" height="541" alt="Pasted image 20261027152459" src="https://github.com/user-attachments/assets/13decefc-775b-41d5-98fd-73ce0daa7bf9" />

This gives one peak at ~0.49 nm and ~3.00 nm. The thickness is the gap between the gaps and the region from ~4 nm to ~6.5 nm, to a total bilayer thickness of ~4.4 nm. Experimental = 4.98 ± 0.15 nm.
___
#### Lateral diffusion

Process the trajectory to remove jumps over the box boundaries with `-pbc nojump`, then calculate the lateral diffusion with `gmx msd`:

```
echo 0 | gmx trjconv -f dspc-eq.xtc -s dspc-eq.tpr -pbc mol -o dspc-eq-pbc.xtc
echo 0 | gmx trjconv -f dspc-eq-pbc.xtc -s dspc-eq.tpr -pbc nojump -o dspc-eq-nojump.xtc
echo 2 | gmx msd -f dspc-eq-nojump.xtc -s dspc-eq.tpr -lateral z -b 15000 -trestart 20 -o msd.xvg
```

<img width="516" height="621" alt="Pasted image 20261027154404" src="https://github.com/user-attachments/assets/7d37a96e-8bd2-42da-bebc-1032ee7a3a44" />

___
#### Order parameters

A script to calculate the (second-rank) order parameter, defined as: $P_2 = 1/2 (3 cos^2<θ> − 1)$, where $θ$ is the angle between the bond and the bilayer normal. $P_2$ = 1 means perfect alignment with the bilayer normal, $P_2$ = −0.5 anti-alignment, and $P_2$ = 0 random orientation, is available [here](https://cgmartini.nl/docs/downloads/tools/legacy.html).

Usage: `do-order-gmx5-p3.py <traj file> <tpr file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <#lipids> <lipid type>`

`python3 do-order-gmx5-p3.py dspc-eq.xtc dspc-eq.tpr 15000 30000 20 0 0 1 128 DSPC`

The script writes the results in `order.dat` 
<img width="785" height="79" alt="Pasted image 20261027160832" src="https://github.com/user-attachments/assets/157357f5-ee2a-4c83-bc29-b968ddb53497" />
and the order parameter profile in `S-profile.dat`, the latter of which we can plot with a custom script:
<img width="633" height="473" alt="Pasted image 20261027161338" src="https://github.com/user-attachments/assets/be445c9e-c640-4301-9774-bfa7197081ea" />
