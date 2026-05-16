# -*- coding: utf-8 -*-
"""
Publication-style RMSD-Rg free-energy landscape plotting.

Input files in --indir:
    FEL_RMSD_Rg_*.csv
Required columns:
    RMSD_nm, Rg_nm, free_energy_kJmol or free_energy_smoothed_kJmol
Optional column:
    count

Outputs:
    plots_unlabeled_publication_style/
    plots_labeled_publication_style/

Default style:
    - each system keeps its own RMSD/Rg axis range
    - free-energy minimum normalized to 0
    - low free energy = deep blue; high free energy = deep red
    - border-connected high-energy/unvisited background forced to vmax
    - contours drawn only on non-background regions
    - both 2D and 3D figures are generated, with and without S1/S2/S3 labels

Example:
    python analysis/plot_fels_publication_style.py --indir ./FELs --sigma 0.35 --vmax 9.31 --bg-frac 0.06
"""

import argparse
import glob
import os
import re
from collections import deque

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def system_name_from_file(path):
    base = os.path.basename(path)
    match = re.match(r"FEL_RMSD_Rg_(.+)\.csv$", base)
    return match.group(1) if match else os.path.splitext(base)[0]


def detect_columns(df, use_mode="energy"):
    lower = {c.lower(): c for c in df.columns}
    rmsd_col = lower.get("rmsd_nm")
    rg_col = lower.get("rg_nm")

    raw_fe = lower.get("free_energy_kjmol")
    smooth_fe = lower.get("free_energy_smoothed_kjmol")
    if use_mode == "smoothed" and smooth_fe is not None:
        fe_col = smooth_fe
    elif raw_fe is not None:
        fe_col = raw_fe
    elif smooth_fe is not None:
        fe_col = smooth_fe
    else:
        fe_col = next((c for c in df.columns if "free_energy" in c.lower()), None)

    count_col = lower.get("count")
    if rmsd_col is None or rg_col is None or fe_col is None:
        raise ValueError(f"Required columns not found. Columns are: {list(df.columns)}")
    return rmsd_col, rg_col, fe_col, count_col


def make_edges(values):
    values = np.asarray(values, dtype=float)
    if values.size == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    mids = 0.5 * (values[:-1] + values[1:])
    first = values[0] - (mids[0] - values[0])
    last = values[-1] + (values[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])


def fel_cmap():
    return LinearSegmentedColormap.from_list(
        "fel_publication",
        [
            "#000080", "#0033cc", "#0066ff", "#00a8ff", "#00d5d5",
            "#56e58e", "#b8ef45", "#ffff33", "#ffbf00", "#ff7700",
            "#ff2300", "#a50000",
        ],
        N=512,
    )


def gaussian_kernel(sigma):
    if sigma <= 0:
        return np.array([1.0])
    radius = max(1, int(3 * sigma))
    x = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-(x * x) / (2 * sigma * sigma))
    return k / k.sum()


def conv_axis(arr, kernel, axis):
    return np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=axis, arr=arr)


def weighted_smooth(z, valid, sigma):
    if sigma <= 0:
        return z.copy()
    k = gaussian_kernel(sigma)
    z0 = np.where(valid, z, 0.0)
    w = valid.astype(float)
    num = conv_axis(conv_axis(z0 * w, k, axis=1), k, axis=0)
    den = conv_axis(conv_axis(w, k, axis=1), k, axis=0)
    return np.divide(num, den, out=np.full_like(num, np.nan), where=(den > 1e-12))


def border_connected(mask):
    mask = np.asarray(mask, dtype=bool)
    ny, nx = mask.shape
    out = np.zeros_like(mask, dtype=bool)
    queue = deque()

    for j in range(nx):
        for i in (0, ny - 1):
            if mask[i, j] and not out[i, j]:
                out[i, j] = True
                queue.append((i, j))
    for i in range(ny):
        for j in (0, nx - 1):
            if mask[i, j] and not out[i, j]:
                out[i, j] = True
                queue.append((i, j))

    while queue:
        i, j = queue.popleft()
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + di, j + dj
            if 0 <= ii < ny and 0 <= jj < nx and mask[ii, jj] and not out[ii, jj]:
                out[ii, jj] = True
                queue.append((ii, jj))
    return out


def load_grid(csv_file, vmax, sigma, bg_frac, use_mode):
    df = pd.read_csv(csv_file)
    rmsd_col, rg_col, fe_col, count_col = detect_columns(df, use_mode)

    grid = df.pivot_table(index=rg_col, columns=rmsd_col, values=fe_col, aggfunc="mean")
    grid = grid.sort_index().sort_index(axis=1)
    x = grid.columns.values.astype(float)
    y = grid.index.values.astype(float)
    z0 = grid.values.astype(float)

    finite = np.isfinite(z0)
    if not np.any(finite):
        raise ValueError(f"No finite free-energy values in {csv_file}")

    z0 = z0 - np.nanmin(z0[finite])
    finite = np.isfinite(z0)
    raw_max = float(np.nanmax(z0[finite]))
    high_cut = raw_max - bg_frac * max(raw_max, 1e-12)
    candidate_bg = (~finite) | (finite & (z0 >= high_cut))

    if count_col is not None:
        try:
            cgrid = df.pivot_table(index=rg_col, columns=rmsd_col, values=count_col, aggfunc="mean")
            cgrid = cgrid.sort_index().sort_index(axis=1)
            count = cgrid.values.astype(float)
            if count.shape == z0.shape:
                candidate_bg = candidate_bg | (np.isfinite(count) & (count <= 0))
        except Exception:
            pass

    bg = border_connected(candidate_bg)
    valid = finite & (~bg)
    zs = weighted_smooth(z0, valid, sigma) if np.any(valid) and sigma > 0 else z0.copy()
    zs = np.where(np.isfinite(zs), zs, z0)
    non_bg = np.isfinite(zs) & (~bg)
    if np.any(non_bg):
        zs = zs - np.nanmin(zs[non_bg])
    else:
        zs = zs - np.nanmin(zs[np.isfinite(zs)])
    zs[bg] = vmax
    zs = np.clip(zs, 0, vmax)
    return x, y, zs, bg, fe_col


def find_basins(x, y, z, bg, n=3, min_sep_x=0.03, min_sep_y=0.015):
    ny, nx = z.shape
    candidates = []
    for i in range(1, ny - 1):
        for j in range(1, nx - 1):
            if bg[i, j]:
                continue
            nb = z[i - 1 : i + 2, j - 1 : j + 2]
            nb_bg = bg[i - 1 : i + 2, j - 1 : j + 2]
            vals = nb[~nb_bg]
            val = z[i, j]
            if vals.size and np.isfinite(val) and val <= np.nanmin(vals):
                candidates.append((float(val), float(x[j]), float(y[i]), j, i))

    if len(candidates) < n:
        candidates = [
            (float(z[i, j]), float(x[j]), float(y[i]), j, i)
            for i in range(ny)
            for j in range(nx)
            if np.isfinite(z[i, j]) and not bg[i, j]
        ]

    candidates.sort(key=lambda t: t[0])
    chosen = []
    for cand in candidates:
        _, xx, yy, _, _ = cand
        if all(abs(xx - old[1]) >= min_sep_x or abs(yy - old[2]) >= min_sep_y for old in chosen):
            chosen.append(cand)
        if len(chosen) >= n:
            break
    return sorted(chosen[:n], key=lambda t: t[1])


def save_basins(outdir, name, basins):
    with open(os.path.join(outdir, f"S1_S2_S3_positions_{name}.csv"), "w", encoding="utf-8") as f:
        f.write("label,RMSD_nm,Rg_nm,Free_energy_kJmol\n")
        for label, item in zip(["S1", "S2", "S3"], basins):
            e, xx, yy, _, _ = item
            f.write(f"{label},{xx:.8f},{yy:.8f},{e:.8f}\n")


def draw_2d(x, y, z, bg, name, outdir, labels_on, vmax):
    cmap = fel_cmap()
    norm = Normalize(vmin=0, vmax=vmax)
    X, Y = np.meshgrid(x, y)
    xe, ye = make_edges(x), make_edges(y)
    z_nonbg = np.ma.array(z, mask=bg)
    levels_filled = np.linspace(0, vmax, 240)
    levels_line = np.linspace(0, vmax, 22)
    basins = find_basins(x, y, z, bg)

    fig, ax = plt.subplots(figsize=(6.8, 5.6), dpi=300)
    ax.set_facecolor("#a50000")
    pm = ax.pcolormesh(xe, ye, z, cmap=cmap, norm=norm, shading="auto")
    ax.contourf(X, Y, z_nonbg, levels=levels_filled, cmap=cmap, norm=norm, extend="max", alpha=0.96)
    ax.contour(X, Y, z_nonbg, levels=levels_line, colors="#d39a00", linewidths=0.55, alpha=0.95)

    if labels_on:
        for label, item in zip(["S1", "S2", "S3"], basins):
            _, xx, yy, _, _ = item
            ax.text(xx, yy, label, color="white", fontsize=18, fontweight="bold", ha="center", va="center")
        save_basins(outdir, name, basins)

    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[0], ye[-1])
    ax.set_xlabel("RMSD (nm)", fontsize=14)
    ax.set_ylabel("Rg (nm)", fontsize=14)
    ax.set_title(f"FEL_RMSD_Rg_{name}", fontsize=16)
    cbar = fig.colorbar(pm, ax=ax, extend="max")
    cbar.set_label("Free energy (kJ/mol)", fontsize=14)
    ticks = np.linspace(0, vmax, 10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.3f}" for t in ticks])

    fig.tight_layout()
    suffix = "labeled" if labels_on else "plain"
    fig.savefig(os.path.join(outdir, f"FEL_RMSD_Rg_{name}_2D_{suffix}.png"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"FEL_RMSD_Rg_{name}_2D_{suffix}.pdf"), bbox_inches="tight")
    plt.close(fig)


def draw_3d(x, y, z, bg, name, outdir, labels_on, vmax):
    cmap = fel_cmap()
    norm = Normalize(vmin=0, vmax=vmax)
    X, Y = np.meshgrid(x, y)
    z_nonbg = np.ma.array(z, mask=bg)
    levels_filled = np.linspace(0, vmax, 240)
    levels_line = np.linspace(0, vmax, 22)
    basins = find_basins(x, y, z, bg)

    fig = plt.figure(figsize=(8.0, 6.3), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, z, rstride=1, cstride=1, cmap=cmap, norm=norm, linewidth=0.10, antialiased=True, shade=False)
    ax.contourf(X, Y, z, zdir="z", offset=0, levels=levels_filled, cmap=cmap, norm=norm)
    ax.contour(X, Y, z_nonbg, zdir="z", offset=0, levels=levels_line, colors="#d39a00", linewidths=0.55)

    if labels_on:
        for label, item in zip(["S1", "S2", "S3"], basins):
            e, xx, yy, _, _ = item
            ax.text(xx, yy, min(e + 0.07 * vmax, 0.92 * vmax), label, fontsize=12, fontweight="bold", color="black")

    ax.set_xlabel("RMSD (nm)", labelpad=8, fontsize=12)
    ax.set_ylabel("Rg (nm)", labelpad=8, fontsize=12)
    ax.set_zlabel("Free energy (kJ/mol)", labelpad=10, fontsize=12)
    ax.set_title(f"FEL_RMSD_Rg_{name}", fontsize=15, pad=12)
    ax.view_init(elev=25, azim=-130)
    ax.set_zlim(0, vmax)
    try:
        ax.set_proj_type("persp")
        xr = float(np.max(x) - np.min(x))
        yr = float(np.max(y) - np.min(y))
        ax.set_box_aspect((xr, yr, 0.70 * vmax))
        ax.zaxis.pane.set_facecolor((0.65, 0.0, 0.0, 1.0))
        ax.xaxis.pane.set_facecolor((1, 1, 1, 1))
        ax.yaxis.pane.set_facecolor((1, 1, 1, 1))
    except Exception:
        pass

    cbar = fig.colorbar(surf, ax=ax, shrink=0.72, pad=0.08, extend="max")
    cbar.set_label("Free energy (kJ/mol)", fontsize=13)
    ticks = np.linspace(0, vmax, 10)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.3f}" for t in ticks])

    fig.tight_layout()
    suffix = "labeled" if labels_on else "plain"
    fig.savefig(os.path.join(outdir, f"FEL_RMSD_Rg_{name}_3D_{suffix}.png"), dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, f"FEL_RMSD_Rg_{name}_3D_{suffix}.pdf"), bbox_inches="tight")
    plt.close(fig)


def process_file(csv_file, out_labeled, out_plain, sigma, vmax, bg_frac, use_mode):
    name = system_name_from_file(csv_file)
    print(f"[INFO] {name}")
    x, y, z, bg, fe_col = load_grid(csv_file, vmax=vmax, sigma=sigma, bg_frac=bg_frac, use_mode=use_mode)
    draw_2d(x, y, z, bg, name, out_plain, labels_on=False, vmax=vmax)
    draw_3d(x, y, z, bg, name, out_plain, labels_on=False, vmax=vmax)
    draw_2d(x, y, z, bg, name, out_labeled, labels_on=True, vmax=vmax)
    draw_3d(x, y, z, bg, name, out_labeled, labels_on=True, vmax=vmax)
    print(f"[DONE] {name}; source column: {fe_col}")


def main():
    parser = argparse.ArgumentParser(description="Publication-style RMSD-Rg FEL plotting.")
    parser.add_argument("--indir", default=".", help="Directory containing FEL_RMSD_Rg_*.csv files")
    parser.add_argument("--sigma", type=float, default=0.35, help="visual smoothing sigma; use 0 for no smoothing")
    parser.add_argument("--vmax", type=float, default=9.31, help="fixed colorbar maximum in kJ/mol")
    parser.add_argument("--bg-frac", type=float, default=0.06, help="border background high-energy threshold fraction")
    parser.add_argument("--use", choices=["energy", "smoothed"], default="energy", help="choose raw or smoothed free-energy column")
    parser.add_argument("--no-smooth", action="store_true", help="disable additional visual smoothing")
    args = parser.parse_args()
    if args.no_smooth:
        args.sigma = 0.0

    indir = os.path.abspath(args.indir)
    os.chdir(indir)
    csv_files = sorted(glob.glob("FEL_RMSD_Rg_*.csv"))
    if not csv_files:
        raise FileNotFoundError("No FEL_RMSD_Rg_*.csv files found.")

    out_labeled = os.path.join(indir, "plots_labeled_publication_style")
    out_plain = os.path.join(indir, "plots_unlabeled_publication_style")
    os.makedirs(out_labeled, exist_ok=True)
    os.makedirs(out_plain, exist_ok=True)

    for csv_file in csv_files:
        process_file(csv_file, out_labeled, out_plain, args.sigma, args.vmax, args.bg_frac, args.use)

    with open(os.path.join(indir, "README_FEL_publication_style.txt"), "w", encoding="utf-8") as f:
        f.write("Publication-style FEL plotting finished.\n")
        f.write(f"input directory: {indir}\n")
        f.write(f"sigma: {args.sigma}\n")
        f.write(f"vmax: {args.vmax}\n")
        f.write(f"bg_frac: {args.bg_frac}\n")
        f.write(f"use: {args.use}\n")
        f.write(f"labeled output: {out_labeled}\n")
        f.write(f"unlabeled output: {out_plain}\n")
        f.write("Free-energy minimum is normalized to 0.\n")
        f.write("Border-connected high-energy/unvisited background is forced to vmax.\n")
        f.write("Contours are drawn only on non-background regions.\n")

    print("\nAll done.")
    print("Labeled output   :", out_labeled)
    print("Unlabeled output :", out_plain)


if __name__ == "__main__":
    main()
