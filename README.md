# Cas12g Molecular Dynamics and FEL Analysis Workflow

This repository contains the final reproducible workflow used for all-atom molecular dynamics simulations and trajectory analysis of Cas12g and selected mutants.

Current target systems:

- Cas12g
- Cas12g F616N
- Cas12g R184S
- Cas12g V623Q
- Cas12g V623S

The repository is intended to store **scripts, parameter templates, methods text and plotting code** only. Large trajectories and generated simulation outputs should not be committed.

## Final workflow overview

The final analysis workflow used in this project is:

```text
protein-only PDB
→ GROMACS topology generation
→ solvation and ion addition
→ energy minimization
→ NVT equilibration
→ Cα-restrained NPT equilibration
→ unrestrained NPT equilibration
→ 500 ns production MD
→ PBC correction and centering
→ first production frame as RMSD reference
→ backbone RMSD / Rg / Cα RMSF / residue-pair distance analysis
→ RMSD-Rg free-energy landscape analysis
→ publication-style 2D and 3D FEL plotting
```

## Simulation settings

The workflow is designed for standard protein-only systems and uses:

- GROMACS
- AMBER99SB-ILDN force field
- TIP3P water model
- Dodecahedral periodic box
- 0.15 M NaCl
- Steepest-descent energy minimization
- 100 ps NVT equilibration
- 1 ns Cα-restrained NPT equilibration
- 500 ps unrestrained NPT equilibration
- 500 ns production MD
- 2 fs integration time step
- LINCS constraints for bonds involving hydrogens
- PME electrostatics
- GPU-accelerated production MD when available

Trajectory analysis and FEL visualization use or resample data at **50 ps = 0.05 ns** intervals. For RMSD, the first frame of the production trajectory is used explicitly as the reference structure.

## Repository contents

```text
Cas12g_md/
├── README.md
├── .gitignore
├── requirements.txt
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
│   ├── plot_fels_publication_style.py
│   ├── plot_rmsd.py
│   └── plot_rmsf.py
└── docs/
    ├── methods_text.md
    └── fel_plotting_guide.md
```

## Quick start for MD

Copy the environment template and edit paths:

```bash
cp env/Cas12g_predict_F616N.env.example Cas12g_predict_F616N.env
vim Cas12g_predict_F616N.env
```

Run a production workflow:

```bash
bash scripts/run_md_new_protein.sh Cas12g_predict_F616N.env
```

For background execution:

```bash
nohup bash scripts/run_md_new_protein.sh Cas12g_predict_F616N.env > pipeline.nohup.out 2>&1 &
echo $! > pipeline.pid
```

Monitor:

```bash
cd /path/to/workdir
watch -n 60 ./md_monitor.sh
tail -f logs/md500.out
```

## Final FEL plotting from existing CSV files

Place the final FEL CSV files in one folder:

```text
FEL_RMSD_Rg_Cas12g.csv
FEL_RMSD_Rg_F616N.csv
FEL_RMSD_Rg_R184S.csv
FEL_RMSD_Rg_V623Q.csv
FEL_RMSD_Rg_V623S.csv
```

Then run:

```bash
python analysis/plot_fels_publication_style.py --indir /path/to/FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.06
```

Outputs:

```text
plots_unlabeled_publication_style/
plots_labeled_publication_style/
README_FEL_publication_style.txt
```

Each system receives separate 2D and 3D FEL plots in PNG and PDF formats. The script also exports S1/S2/S3 basin coordinates for the labeled plots.

## Required Python packages

```bash
pip install -r requirements.txt
```

Required packages:

- numpy
- pandas
- matplotlib

## Important notes

This workflow is for standard amino-acid protein-only systems. If the input contains ligands, cofactors, metals, nucleic acids, glycans, nonstandard residues or membrane components, additional topology preparation is required.

Do **not** commit large generated files such as:

- `*.xtc`, `*.trr`, `*.tng`
- `*.edr`, `*.cpt`
- production `*.gro` or `*.tpr`
- large logs or nohup outputs
- complete MD working directories

Use external storage such as Zenodo, Figshare, OSF, institutional storage or data-on-request for production trajectories.

## Citation-ready method summary

All-atom molecular dynamics simulations were performed using GROMACS with the AMBER99SB-ILDN force field and TIP3P water model. Protein-only structures were placed in dodecahedral periodic boxes, solvated, neutralized and supplemented with 0.15 M NaCl. Systems were energy-minimized using the steepest-descent algorithm, followed by 100 ps NVT equilibration, 1 ns Cα-restrained NPT equilibration and 500 ps unrestrained NPT equilibration. Production simulations were designed for 500 ns with a 2 fs integration time step. Trajectories were corrected for periodic boundary artifacts, centered and fitted using backbone atoms. Backbone RMSD was calculated relative to the first production frame. Rg, Cα RMSF and selected residue-pair distances were calculated from the processed trajectories. RMSD-Rg free-energy landscapes were constructed from 50 ps interval data using relative free energy, with the minimum normalized to 0 kJ/mol. Unvisited or low-probability regions were assigned to the highest free-energy level and displayed in deep red.
