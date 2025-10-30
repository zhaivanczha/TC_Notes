## Step 1: Generate orbital (molden $\rightarrow$ gwfn.data)
#### package: 
- PySCF / Molpro
#### installation: 
- pip install pyscf / module load molpro/2015.1 (on `alamdlogin1/2`)
#### usage: 
- Run jobs as documented in PySCF/Molpro (One can find `create_molden_file.py` for PySCF in the utils directory in `PYTCHINT` package)
- Note: 
	- the molden will be changed to `gpcc.casl` in the `create_molden_file.py` using the script `$CASINO/../utils/wfn_converters/molden/molden2qmc.py`, so **CASINO** have to be pre-installed
	- a FCIDUMP is also needed for Step 3 in PYTCHINT in case `use_fcidump true` 
#### important parameter settings: 
1. Basis set
2. Orbital type
	- HF
	- DFT
	- MP2
	- SA-CASSCF(multi-determinant optimization)
3. Open-shell
	- RHF/ROHF/UHF
4. Point group symmetry
## Step 2: Generate Jastrow
#### package: 
- CASINO[TODO: deterministic optimization in PYTCHINT]
#### installation: 
1. git clone git@github.com:fkfest/CASINO.git (private repository of FKFEST)
2. cd CASINO && ./install
3. choose `git@github.com:fkfest/CASINO.git` and choose the one with `matches hostname explicitly`
4. choose `Compile CASINO for already-configured CASINO_ARCHs`  and type the number of CASINO_ARCH you would like to compile (you can choose no feature)
5. choose `Save (unmodified) configuration to .bashrc.casino file and quit`
- Note: if using CASINO on PC, one should have an MPI-enabled gfortran compiler in your path and choose the `linuxpc-gcc-parallel` instead of the recommended `linuxpc-gcc-parallel.fkf` (outdated)
#### usage:
The usage of CASINO is a bit difficult, it is highly recommended that one use the scripts `01_gen_gpcc.sh`,`02a_opt.sh`,... with python scripts (see `scripts` folder; adopted from @Philip Haupt) to automatically generate the input files. Though, one can write the input file like this:
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
and change it accordingly. Then run the command `runqmc -n1 -c -xxx` to generate the `casino` slurm submitting script (the `-c` means only generate script without submission).
In principle, there are two `runtype` calculations in Jastrow optimization of CASINO for TC: 
1. `gen_gpcc/gen_gpcc_simple/gen_gpcc_simple` (~~optional~~/mandatory)  [TODO: difference between them; if this will influence the en cusp in $\chi$ ]
	- A Jastrow factor to correct the electron-nucleui cusp (e-n cusp) for orbitals: 
		   $f(r) = \left[ \log \left( e^{\sum_{k=0}^{4} \alpha_k r^k} + C \right) - \log \phi_s(r) \right] \Theta(L - r)$
	- Also called `hat` function; it is used to stabilize the VMC optimization procedure.
	- Before the next step, the generated `gpcc.gasl` should be copied to the `vmc_opt` directory and renamed as `parameters.casl`
	- It seems that now the CASINO will default assume the user use the gpcc and the determined parameter in $\chi$ will no longer consider the cusp condition (refer to [[TC_workflow#^1f73cb]]). 
2. `vmc_opt` (mandatory)
	- There should be two steps of VMC, the first is a rough VMC followed by a more precise one
	- The `gjastrow` form of CASINO can be used to generate different types of Jastrow factor. The commonly-used one is Drummond–Towler–Needs (DTN) factor, which can be expanded by `natural power` basis function: $\phi_k(r)=r^{k-1}$
     times `polynomial` cutoff function: $f(r)=(1-r/L)^C\Theta(L-r)$.
     The Jastrow factor is written as:  
   $$J(r)=e^{[\sum_k c_k \phi_k(r)]\times f(r)}$$.  
     Therefore, the e-e cusp condition $\lim_{r_{ij}\rightarrow 0} \frac{\partial u(r_{ij})}{\partial r_{ij}}=1/2$ will result in:
   $$c_1 \times(-\frac{C}{L})+ c_2=1/2 \rightarrow c_1=(c_2-1/2)\times\frac{L}{C}$$  
   ; the same, the e-n cusp condition $\lim_{r_{Ii}\rightarrow 0} \frac{\partial \chi(r_{Ii})}{\partial r_{Ii}}=-Z$ result in:
   $$c_1=(c_2+Z)\times\frac{L}{C}.$$  
     However, consideration of gpcc has already treated the en cusp, as a result, $\lim_{r_{Ii}\rightarrow 0} \frac{\partial \chi(r_{Ii})}{\partial r_{Ii}}=0$ and $c_1=c_2\times\frac{L}{C}$.  
     The CASINO will not print the determined $c_1$ out unless add `Print determined: T` in the `parameters.casl`. ^1f73cb
	- There are different channels of `gjastrow`, it can be used to distinguish the different spin of electron pairs; the e-n terms of different atoms. Usually, we only distinguish the e-n terms between different elements.
	- `opt_method` should be set to `varmin` to avoid the over-corrected cusp effect, the `emin` will introduce the non-variational energy in the TC result. (`emin` in VMC can provide a lower VMC energy and therefore DMC energy, *but not stable as `varmin`*[TODO: not sure])
- Note:
	1. The generally used parameters for DTN Jastrow is:
		- e-e order: 9, e-e cutoff: 4.5
		- e-n order: 9, e-n cutoff: 4 (1)
		- e-e-n order: 4; e-e-n cutoff: 4 (2)
	2. The Jastrow factor used in VMC and TC only rely on the radial distance between two particles (isotropic), i.e. the r is a scalar instead of vector. (not important and also brings a lot computational efforts)
	3. `DMC` will be contained in all files, even though it will not be performed
	4. Backflow should not be used in Jastrow optimization for TC, since it gives optimal Jastorw factor for backflow Slater-jastrow but not for Slater-Jastrow
	5. The speed and memory requirement of VMC is not related to the basis set, so small basis set might have a slower speed and also larger memory requirement. E.g. STO-3G for $N_2$ and cc-pvdz for $CO_2$
	6. There are some problems with the nodes on `alamdlogin2`: 
		- `gen2` needs adding `export PMIX_MCA_psec=^munge`;
		- `gen2.big` needs changing `module load gnu-openmp` to `module load software/gnu-openmpi/4.1.5` in the input file `casino`
	 7. TODO: it is possible to set `opt_det_coeff = T` for TC, a better reference?
#### important parameter settings: 
1. Jastrow form
	- Type:
		- DTN
			- Order
			- Cutoff
		- Boys–Handy
		- Composed form by hand
		- ...
	- Different channel choice
		- 1=2,Z
		- ...
2. Optimization way
	- varmin
	- emin
	- madmin
3. Optimization step (can be determined automatically by the scripts)
## Step 3: Integral evaluation
#### package:
PYTCHINT
#### installation:
1. git clone git@github.com:fkfest/pytchint.git (private repository of FKFEST) && cd pytchint
2. mamba env create -f environment.yml  # or: conda env create -f environment.yml
3. conda activate pytchint
4. pip install -e . (editable)
5. python3 test_runner.py # pytest (now some cases fails)
#### usage:
One can generate a input file like this:
``` tchint.yml
nelec: 7
spin: 1
grid_lvl: 2
eval_mode_str: xtc
xtc_perf_level: 1
memory_logfile: memlog
jastrow_id: casino
grid_blocking: false

fcidump_filename_in: FCIDUMP.bare
use_fcidump: true


combined_jastrow: true
verbose: 5
```
here the `FCIDUMP.bare` is the generated `FCIDUMP` from Step 1. And then one can use the `tchint_standalone.py` in the utils directory in `PYTCHINT` package to run PYTCHINT: `mpirun -np {ntasks} ../tchint_standalone.py --debug-log tchint.yml`.

#### important parameter settings: 
1. eval_mode_str
	- xTC
	- full TC
2. grid_lvl [TODO: is this for the evaluation of K/L tensors?]
## Step 4: quantum chemistry methods
#### package:
ElemCo.jl / NECI / Block2/ ...
#### installation:
refer to the document pls. x_x
