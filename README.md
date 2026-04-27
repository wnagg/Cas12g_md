# Cas12g Molecular Dynamics Workflow

This repository contains a reproducible GROMACS workflow for all-atom molecular dynamics simulations of the Cas12g F616N mutant protein.

The workflow is designed for a **protein-only** system prepared from `Cas12g.pdb`, using:

- GROMACS 2024.5
- AMBER99SB-ILDN force field
- TIP3P water model
- Dodecahedral periodic box
- 0.15 M NaCl
- Steepest-descent energy minimization
- C-alpha restrained NVT equilibration
- C-alpha restrained NPT equilibration
- Unrestrained NPT relaxation
- GPU-accelerated production MD
- Lightweight trajectory output for long simulations

## Repository contents

```text
Cas12g_md/
├── README.md
├── .gitignore
├── env/
│   └── Cas12g_predict_F616N.env.example
├── scripts/
│   └── run_md_new_protein.sh
├── mdp/
│   ├── ions.mdp
│   ├── minim.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   ├── npt_500ps_noca.mdp
│   └── md_500ns.mdp
├── analysis/
│   ├── plot_rmsd.py
│   └── plot_rmsf.py
└── docs/
    ├── operation_guide.md
    └── methods_text.md
```

## Important notes

This workflow is intended for standard amino-acid protein systems. If the input structure contains ligands, metal ions, cofactors, modified residues, nucleic acids, glycans, or membrane components, additional topology preparation is required before running this workflow.

Large files should not be committed to GitHub. Do **not** upload production trajectories such as `md500.xtc`, `md500.trr`, full solvated coordinate files, checkpoint files, or large energy files. Use Zenodo, Figshare, OSF, institutional storage, or external data availability on request for large simulation data.

## Quick start

Copy your protein-only PDB to the server, then edit the example environment file:

```bash
cp env/Cas12g_predict_F616N.env.example Cas12g_predict_F616N.env
vim Cas12g_predict_F616N.env
```

For a short test run, set:

```bash
MD_NS=1
XTC_INTERVAL_PS=10
RUN_LIGHT_ANALYSIS=no
```

Start the workflow:

```bash
chmod +x scripts/run_md_new_protein.sh
nohup bash scripts/run_md_new_protein.sh Cas12g_predict_F616N.env > pipeline_1ns.nohup.out 2>&1 &
echo $! > pipeline_1ns.pid
```

Monitor:

```bash
cd /path/to/workdir
watch -n 60 ./md_monitor.sh
tail -f logs/md500.out
```

After the 1 ns test completes without fatal errors, LINCS warnings, or NaN values, prepare a 500 ns production configuration:

```bash
MD_NS=500
XTC_INTERVAL_PS=100
RUN_LIGHT_ANALYSIS=yes
ANALYSIS_BEGIN_NS=50
```

## Citation-ready method summary

All-atom MD simulations were performed using GROMACS 2024.5 with the AMBER99SB-ILDN force field and TIP3P water model. The protein was placed in a dodecahedral periodic box with a minimum solute-box distance of 1.0 nm, solvated, neutralized, and supplemented to 0.15 M NaCl. The system was energy-minimized with the steepest-descent algorithm, followed by C-alpha restrained NVT and NPT equilibration and an additional unrestrained NPT relaxation. Production MD was performed under the NPT ensemble with a 2 fs time step, LINCS constraints on bonds involving hydrogens, PME electrostatics, V-rescale temperature coupling, and C-rescale pressure coupling. Compressed coordinates were saved at user-defined intervals to reduce storage requirements.
