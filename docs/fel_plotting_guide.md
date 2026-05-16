# FEL plotting guide

This guide documents the final publication-style FEL plotting approach used in this project.

## Input

Place final FEL CSV files in one directory:

```text
FEL_RMSD_Rg_Cas12g.csv
FEL_RMSD_Rg_F616N.csv
FEL_RMSD_Rg_R184S.csv
FEL_RMSD_Rg_V623Q.csv
FEL_RMSD_Rg_V623S.csv
```

Each CSV should contain at least:

```text
RMSD_nm,Rg_nm,free_energy_kJmol
```

The script also accepts `free_energy_smoothed_kJmol` and optional `count`.

## Run

From the repository root:

```bash
python analysis/plot_fels_publication_style.py --indir /path/to/FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.06
```

For Windows, use the same command in PowerShell or Anaconda Prompt:

```powershell
python analysis\plot_fels_publication_style.py --indir D:\FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.06
```

## Output

The script creates:

```text
plots_unlabeled_publication_style/
plots_labeled_publication_style/
README_FEL_publication_style.txt
```

Each system receives:

```text
FEL_RMSD_Rg_<system>_2D_plain.png
FEL_RMSD_Rg_<system>_2D_plain.pdf
FEL_RMSD_Rg_<system>_3D_plain.png
FEL_RMSD_Rg_<system>_3D_plain.pdf
FEL_RMSD_Rg_<system>_2D_labeled.png
FEL_RMSD_Rg_<system>_2D_labeled.pdf
FEL_RMSD_Rg_<system>_3D_labeled.png
FEL_RMSD_Rg_<system>_3D_labeled.pdf
S1_S2_S3_positions_<system>.csv
```

## Parameters

- `--sigma`: visual smoothing intensity. Default: `0.35`. Use `0.0` or `--no-smooth` for no extra visual smoothing.
- `--vmax`: fixed maximum of the color scale in kJ/mol. Default: `9.31`.
- `--bg-frac`: fraction used to identify high-energy border-connected background. Default: `0.06`.
- `--use`: choose `energy` or `smoothed`. Default: `energy`.

Recommended initial command:

```bash
python analysis/plot_fels_publication_style.py --indir /path/to/FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.06
```

If the background is not red enough:

```bash
python analysis/plot_fels_publication_style.py --indir /path/to/FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.10
```

If the plot looks over-smoothed:

```bash
python analysis/plot_fels_publication_style.py --indir /path/to/FELs --sigma 0.2 --vmax 9.31 --bg-frac 0.06
```

## Methodological note

The plotted free energy is relative free energy. The minimum free-energy state is normalized to 0 kJ/mol. Border-connected high-energy or unvisited regions are assigned to the fixed color-scale maximum and therefore appear as a deep-red high-energy background. Contours are drawn only over non-background regions to avoid artificial yellow edges around the plot.
