# Methods text for manuscript

## Molecular dynamics simulations

All-atom molecular dynamics simulations were performed for Cas12g WT and the F616N, R184S, V623Q and V623S single-point mutants using GROMACS. Each mutant was generated independently from the same WT PDB by PyMOL side-chain mutagenesis; the WT coordinates were not overwritten. Mutation identities and residue numbering were checked before topology generation. Protein-only structures were parameterized with the AMBER99SB-ILDN force field and TIP3P water model. Each system was placed in a dodecahedral periodic box, neutralized and supplemented with NaCl to 0.15 M.

Energy minimization used steepest descent. Equilibration comprised 100 ps NVT, 1 ns NPT with C-alpha positional restraints and 500 ps unrestrained NPT. Temperature and pressure were maintained at 300 K and 1 bar. Bonds involving hydrogen atoms were constrained using LINCS, and long-range electrostatics were treated with particle-mesh Ewald. Production MD used a 2 fs timestep and GPU acceleration when supported. The default panel configuration uses 100 ns per system for exploratory screening; replace this with the actual duration of the trajectories used in the final analysis. Short screening runs alone should not be described as demonstrating equilibrium convergence.

## Trajectory preprocessing and RMSD, Rg and RMSF analyses

Periodic-boundary artifacts were removed, proteins were centered, and frames were least-squares fitted using backbone atoms. The first frame of each production trajectory was used as the RMSD reference. Backbone RMSD and radius of gyration (Rg) were calculated throughout the trajectory. Per-residue C-alpha RMSF was calculated after excluding the early production segment; the default analysis starts at 20 ns. The default trajectory output interval was 50 ps. Report the actual interval and analysis window used.

## Residue-pair distance analysis

Distances between R661-R224, D513-E655 and D513-K728 were calculated from PBC-corrected, centered and backbone-fitted trajectories. The implemented selections measure the R661 guanidinium group to R224 CZ, the D513 and E655 carboxylate groups, and the D513 carboxylate group to K728 NZ. Confirm residue numbering against the processed topology if the input structure is renumbered or has missing residues.

## RMSD-Rg free-energy landscape analysis

The two-dimensional probability distribution of backbone RMSD and Rg was converted to relative free energy using Delta G = -kBT ln P, where P is the probability of an occupied histogram bin, kB is the Boltzmann constant and T is 300 K. The minimum observed free energy was normalized to zero. Empty bins were masked as unobserved rather than treated as measured free-energy values. The comparison panel uses common RMSD/Rg axes and a shared color scale. These landscapes are descriptive estimates; their apparent basins and relative depths should not be interpreted as converged thermodynamics without adequate sampling and convergence checks.

## Initial structure quality control

Because mutants are constructed from WT coordinates by side-chain mutagenesis, initial global C-alpha/backbone RMSDs after alignment are expected to be nearly zero. These values verify residue numbering and preservation of the starting backbone; they are not evidence of mutation-induced backbone rearrangement.

## Machine-specific reporting

Record the exact GROMACS version, CPU and GPU using gmx --version, nvidia-smi and lscpu. Report the actual production duration, sampling interval, analysis window and number of independent replicas. Do not describe a 100-ns screening run as a 500-ns production run.
