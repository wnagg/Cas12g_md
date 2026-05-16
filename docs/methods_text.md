# Methods text for manuscript

## Molecular dynamics simulations

All-atom molecular dynamics simulations were performed for Cas12g and selected mutants, including F616N, R184S, V623Q and V623S, using GROMACS. The initial structures were protein-only PDB files with water molecules, exogenous ligands and hydrogen atoms removed. Protein topologies were generated using the AMBER99SB-ILDN force field, and the TIP3P water model was used for solvation. Each system was placed in a dodecahedral periodic water box, neutralized with counterions and supplemented with NaCl to a final concentration of 0.15 M.

Energy minimization was performed using the steepest-descent algorithm until the maximum force was below 1000 kJ mol^-1 nm^-1. The systems were then equilibrated by 100 ps NVT equilibration, followed by 1 ns NPT equilibration with C-alpha positional restraints and 500 ps unrestrained NPT equilibration. Production simulations were designed for 500 ns with a 2 fs integration time step. The temperature and pressure were maintained at 300 K and 1 bar, respectively. Bonds involving hydrogen atoms were constrained using the LINCS algorithm, and long-range electrostatic interactions were treated using the particle mesh Ewald method. GPU acceleration was used for production simulations when available.

## Trajectory preprocessing

Prior to trajectory analysis, periodic boundary artifacts were removed and proteins were centered. The first frame of each production trajectory was extracted as the reference structure. All frames were least-squares fitted based on backbone atoms before RMSD-related analyses. Therefore, backbone RMSD values represent conformational deviations relative to the initial production conformation.

## RMSD, Rg and RMSF analyses

Backbone RMSD was calculated relative to the first frame of the production trajectory. The radius of gyration was calculated to evaluate protein compactness. Per-residue C-alpha RMSF was calculated using trajectory frames after 50 ns of production simulation to reduce the influence of early equilibration. Trajectory data were analyzed or resampled at a 50 ps interval for comparative visualization and RMSD-Rg free-energy landscape construction.

## Residue-pair distance analysis

Distances between selected residue pairs were monitored throughout the production trajectories to evaluate mutation-associated changes in local interaction networks. The analyzed residue pairs included R661-R224, D513-E655 and D513-K728. Distances were extracted from PBC-corrected, centered and backbone-fitted trajectories.

## Free-energy landscape analysis

Free-energy landscapes were constructed using backbone RMSD and Rg as reaction coordinates. The two-dimensional probability distribution of RMSD and Rg was converted into relative free energy according to:

```text
Delta G = -kBT ln P
```

where P is the probability of each conformational state, kB is the Boltzmann constant and T was set to 300 K. The minimum free-energy state was normalized to 0 kJ/mol. Unvisited or low-probability regions were assigned to the highest free-energy level and displayed in deep red. Publication-style two-dimensional and three-dimensional FEL plots were generated using Python with NumPy, Pandas and Matplotlib. Major low-energy basins were identified from local minima on the RMSD-Rg FEL maps and labeled as S1, S2 and S3.

## Recommended wording

Production simulations were designed for 500 ns with a 2 fs integration time step. Trajectory frames were analyzed or resampled at a 50 ps interval for RMSD-Rg free-energy landscape construction. The first production frame was used as the reference structure for backbone RMSD calculation.

Before submission, verify and insert the exact GROMACS version, CPU model and GPU model using:

```bash
gmx --version
nvidia-smi
lscpu
```
