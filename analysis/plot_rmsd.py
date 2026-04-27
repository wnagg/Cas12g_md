#!/usr/bin/env python3
"""Plot backbone RMSD from a two-column GROMACS XVG-derived data file or CSV."""
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

infile = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('rmsd_backbone.csv')
outfile = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('rmsd_backbone.png')

try:
    data = np.loadtxt(infile, delimiter=',', skiprows=1)
except Exception:
    data = np.loadtxt(infile)

t = data[:, 0]
rmsd = data[:, 1]
plt.figure(figsize=(8, 5))
plt.plot(t, rmsd, lw=0.8)
plt.xlabel('Time (ns)')
plt.ylabel('Backbone RMSD (nm)')
plt.tight_layout()
plt.savefig(outfile, dpi=300)
print(f'Wrote {outfile}')
