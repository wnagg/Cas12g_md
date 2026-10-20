# Cas12g Molecular Dynamics and FEL Analysis Workflow

This repository now includes a one-command panel workflow for Cas12g WT and four single-point mutants: F616N, R184S, V623Q and V623S. Each mutant is independently constructed from the original WT PDB using PyMOL, followed by separate GROMACS simulations and trajectory analysis.

## Quick start

Copy env/Cas12g_panel.env.example to env/Cas12g_panel.env. Edit GMX, PYMOL, WT_PDB, CHAIN and WORKDIR. Use the original WT PDB, not a mutant.

Run:
    bash scripts/run_cas12g_panel.sh env/Cas12g_panel.env

Background execution:
    nohup bash scripts/run_cas12g_panel.sh env/Cas12g_panel.env > cas12g_panel.nohup.out 2>&1 &
    echo $! > cas12g_panel.pid

Modes:
    bash scripts/run_cas12g_panel.sh env/Cas12g_panel.env --build-only
    bash scripts/run_cas12g_panel.sh env/Cas12g_panel.env --run-only
    bash scripts/run_cas12g_panel.sh env/Cas12g_panel.env --analyze-only

The five systems run sequentially to avoid competing for GPU memory. Existing structure files are reused. Set REBUILD_MUTANTS=yes only when intentionally rebuilding; do not rebuild structures in a work directory containing production TPR files.

## Resource-conscious defaults

- Integration timestep: 2 fs with hydrogen-bond constraints; no unvalidated 4-fs timestep.
- NVT: 100 ps; C-alpha-restrained NPT: 1 ns; unrestrained NPT: 500 ps.
- Production: 100 ns per system for exploratory screening.
- Trajectory output: every 50 ps, about 2,000 frames per 100 ns.
- RMSF/FEL analysis begins at 20 ns; FEL bin count is limited for short trajectories.

A 100-ns trajectory is a screening run, not proof of equilibrium convergence or converged free-energy landscapes. Extend relevant systems to 200–500 ns or longer and consider independent replicas for final mechanistic claims.

## Outputs

WORKDIR contains structures/ (WT and mutants), runs/<variant>/ (separate GROMACS run and analysis directories), and analysis/ (initial_models_qc.csv, comparison_timeseries.pdf/png, comparison_FEL.pdf/png). Per-system analysis includes backbone RMSD relative to the first production frame, Rg, C-alpha RMSF, three selected residue-pair distances, an RMSD-Rg FEL, CSV tables and a summary figure. Processed intermediate trajectories and XVG files are temporary and removed automatically.

The initial-model QC checks mutation identity and fitted C-alpha/backbone RMSD. Since mutations are made by side-chain mutagenesis on WT coordinates, near-zero backbone RMSD is expected; this is construction QC, not evidence of a mutation-induced global structural change. Empty FEL bins are masked as unobserved.

Monitor progress with:
    tail -f /path/to/WORKDIR/runs/WT/logs/10_md500.out
    watch -n 60 /path/to/WORKDIR/runs/WT/md_monitor.sh

Rerun the same command and config to resume incomplete stages/checkpoints. Do not change production duration, timestep or starting structure inside an existing work directory.

Python dependencies: numpy, pandas and matplotlib. Install with pip install -r requirements.txt. Record GROMACS/CPU/GPU details using gmx --version, nvidia-smi and lscpu. Report the actual production duration and analysis window in the manuscript; do not report a planned 500-ns protocol if only 100 ns was analyzed.
