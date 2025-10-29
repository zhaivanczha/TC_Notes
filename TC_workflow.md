## Step 1: Generate orbital (molden $\rightarrow$ gwfn.data)
#### package: 
- PySCF / Molpro
#### installation: 
- pip install pyscf / module load molpro/2015.1 (on `alamdlogin1/2`)
#### usage: 
- run jobs as documented in PySCF/Molpro (One can find `create_molden_file.py` for PySCF in the utils directory in `PYTCHINT` package)
- Note that: the molden will be changed to `gpcc.casl` in the `create_molden_file.py` using the script `$CASINO/../utils/wfn_converters/molden/molden2qmc.py`, so **CASINO** have to be pre-installed
#### parameter settings: 
1. Basis set
2. Orbital type
	- HF
	- DFT
	- MP2
	- SA-CASSCF(multi-determinant optimization)
3. Open-shell
	- RHF/ROHF/UHF
4. point group symmetry
## Step 2: Generate Jastrow
#### package: 
- CASINO
#### installation: 
1. git clone git@github.com:fkfest/CASINO.git (private repository of FKFEST)
2. cd CASINO && ./install
3. choose `git@github.com:fkfest/CASINO.git` and choose the one with `matches hostname explicitly`
4. choose `Compile CASINO for already-configured CASINO_ARCHs`  and type the number of CASINO_ARCH you would like to compile (you can choose no feature)
5. choose `Save (unmodified) configuration to .bashrc.casino file and quit`
- Note: if using CASINO on PC, one should have an MPI-enabled gfortran compiler in your path and choose the `linuxpc-gcc-parallel` instead of the recommended `linuxpc-gcc-parallel.fkf` (outdated)
#### usage:
The usage of CASINO is a bit difficult, one can write the input file like this:
```
#-------------------#
# CASINO input file #
#-------------------#

# SYSTEM
neu               : $neu           #*! Number of up electrons (Integer)
ned               : $ned           #*! Number of down electrons (Integer)
periodic          : F              #*! Periodic boundary conditions (Boolean)
atom_basis_type   : gaussian       #*! Basis set type (text)
psi_s             : slater         #*! Type of [anti]symmetrizing wfn (Text)
complex_wf        : F              #*! Wave function real or complex (Boolean)
cusp_correction   : F
use_gpcc          : F

# RUN
runtype           : vmc_opt        #*! Type of calculation (Text)
newrun            : T              #*! New run or continue old (Boolean)
testrun           : F              #*! Test run flag (Boolean)

# VMC
vmc_equil_nstep   : 5000           #*! Number of equilibration steps (Integer)
vmc_nstep         : 1000000        #*! Number of steps (Integer)
vmc_nblock        : 10             #*! Number of checkpoints (Integer)
vmc_nconfig_write : 100000         #*! Number of configs to write (Integer)
vmc_sample_hf     : T
vmc_decorr_period : 50

# DMC
dmc_equil_nstep   : 2000           #*! Number of steps (Integer)
dmc_equil_nblock  : 1              #*! Number of checkpoints (Integer)
dmc_stats_nstep   : 10000          #*! Number of steps (Integer)
dmc_stats_nblock  : 1              #*! Number of checkpoints (Integer)
dmc_target_weight : 1000.0         #*! Total target weight in DMC (Real)
dtdmc             : 0.002          #*! DMC time step (Real)

# OPTIMIZATION
opt_method        : varmin         #*! Opt method (varmin/madmin/emin/...)
opt_cycles        : 8              #*! Number of optimization cycles (Integer)
opt_jastrow       : T              #*! Optimize Jastrow factor (Boolean)
opt_det_coeff     : F              #*! Optimize determinant coeffs (Boolean)
opt_backflow      : F              #*! Optimize backflow parameters (Boolean)

%block opt_plan
1 maxiter=1 method=madmin sample_hf=F
2 maxiter=1 method=madmin
3 maxiter=1
4 maxiter=1
5
6
7
8
%endblock opt_plan

# GENERAL PARAMETERS
use_jastrow       : T              #*! Use a Jastrow function (Boolean)
expot             : $expot
backflow          : F              #*! Use backflow corrections (Boolean)
timing_info       : F              #*! Activate subroutine timers (Boolean)
checkpoint : -1
```
, and change it accordingly. Then run the command `runqmc -n1 -c -xxx` to generate the `casino` slurm submitting script (the `-c` means only generate script without submission).
In principle, there are two `runtype` calculations in Jastrow optimization of CASINO: 
1. `gen_gpcc/gen_gpcc_simple/gen_gpcc_simple` (optional)  [TODO: difference between them; the reason why using them]
	- A Jastrow factor to correct the electron-nucleui cusp (e-n cusp) for orbitals: 
		   $f(r) = \left[ \log \left( e^{\sum_{k=0}^{4} \alpha_k r^k} + C \right) - \log \phi_s(r) \right] \Theta(L - r)$
	- also called `hat` function
	- Before the next step, the generated `gpcc.gasl` should be copied to the `vmc_opt` directory and renamed as `parameters.casl`
2. `vmc_opt` (mandatory)
	- there should be two steps of VMC, the first is a rough VMC followed by a more precise one
	- 
- Note:
	1. The generally used parameters for DTN Jastrow is:
		- e-e order: 9, e-e cutoff: 4.5
		- e-n order: 9, e-n cutoff: 4 (1)
		- e-e-n order: 4; e-e-n cutoff: 4 (2)
	2. `DMC` will be contained in all files, even though it will not be performed
	3. Backflow should not be used in Jastrow optimization for TC, since it gives optimal Jastorw factor for backflow Slater-jastrow but not for Slater-Jastrow
	4. The speed and memory requirement of VMC is not related to the basis set, so small basis set might have a slower speed and also larger memory requirement. E.g. STO-3G for $N_2$ and cc-pvdz for $CO_2$
	5. There are some problems with the nodes on `alamdlogin2`: 
		- `gen2` needs adding `export PMIX_MCA_psec=^munge`;
		- `gen2.big` needs changing `module load gnu-openmp` to `module load software/gnu-openmpi/4.1.5` in the input file `casino`
#### settings:
1. 