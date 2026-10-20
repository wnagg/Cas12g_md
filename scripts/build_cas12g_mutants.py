#!/usr/bin/env python3
"""Build Cas12g single-point mutants from the same WT PDB using PyMOL mutagenesis.

Run inside PyMOL, for example:
  CAS12G_WT_PDB=/path/to/Cas12g_predict.pdb \
  CAS12G_STRUCT_DIR=/path/to/workdir/structures \
  pymol -cq -d 'run scripts/build_cas12g_mutants.py' -d 'quit'

All variants are independently loaded from WT. No mutant is used as the template
for another mutant. The WT PDB is never modified in place.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

from pymol import cmd

# tag: (residue number, WT residue name, mutant residue name)
MUTATIONS = {
    "F616N": (616, "PHE", "ASN"),
    "R184S": (184, "ARG", "SER"),
    "V623Q": (623, "VAL", "GLN"),
    "V623S": (623, "VAL", "SER"),
}


def residue_at(pdb_path: Path, chain: str, resi: int) -> str | None:
    """Read the residue identity from an ATOM CA record, independent of PyMOL."""
    found: set[str] = set()
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        if line[12:16].strip() != "CA":
            continue
        if (line[21].strip() or "_") != (chain.strip() or "_"):
            continue
        if line[22:26].strip() != str(resi) or line[26].strip():
            continue
        found.add(line[17:20].strip().upper())
    if len(found) > 1:
        raise RuntimeError(f"Ambiguous CA residue identity at chain {chain}, residue {resi}: {found}")
    return next(iter(found)) if found else None


def main() -> None:
    wt_env = os.environ.get("CAS12G_WT_PDB", "").strip()
    out_env = os.environ.get("CAS12G_STRUCT_DIR", "").strip()
    chain = os.environ.get("CAS12G_CHAIN", "A").strip() or "A"
    if not wt_env or not out_env:
        raise RuntimeError("Set CAS12G_WT_PDB and CAS12G_STRUCT_DIR before launching PyMOL.")

    wt_path = Path(wt_env).expanduser().resolve()
    out_dir = Path(out_env).expanduser().resolve()
    if not wt_path.is_file():
        raise FileNotFoundError(f"WT PDB not found: {wt_path}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Validate all WT identities before making any changes.
    for tag, (resi, expected, target) in MUTATIONS.items():
        observed = residue_at(wt_path, chain, resi)
        if observed != expected:
            raise RuntimeError(
                f"Mutation preflight failed for {tag}: expected WT {expected} at "
                f"chain {chain}, residue {resi}, found {observed!r}. Check PDB numbering/chain."
            )

    shutil.copy2(wt_path, out_dir / "WT.pdb")
    print(f"[OK] WT copied unchanged: {out_dir / 'WT.pdb'}")

    cmd.feedback("disable", "all", "actions")
    for tag, (resi, expected, target) in MUTATIONS.items():
        out_path = out_dir / f"{tag}.pdb"
        cmd.reinitialize()
        cmd.load(str(wt_path), "wt")

        selection = f"wt and chain {chain} and resi {resi} and name CA"
        model = cmd.get_model(selection)
        if len(model.atom) != 1:
            raise RuntimeError(
                f"PyMOL found {len(model.atom)} CA atoms for chain {chain}, residue {resi}; expected exactly one."
            )
        observed = model.atom[0].resn.upper()
        if observed != expected:
            raise RuntimeError(f"PyMOL identity mismatch for {tag}: expected {expected}, found {observed}")

        cmd.wizard("mutagenesis")
        wizard = cmd.get_wizard()
        wizard.set_mode(target)
        wizard.do_select(f"wt and chain {chain} and resi {resi}")
        wizard.apply()
        cmd.set_wizard()
        cmd.save(str(out_path), "wt", state=1)
        cmd.delete("all")

        output_residue = residue_at(out_path, chain, resi)
        if output_residue != target:
            raise RuntimeError(
                f"PyMOL output validation failed for {tag}: expected {target} at "
                f"chain {chain}, residue {resi}; found {output_residue!r}."
            )
        if not out_path.is_file() or out_path.stat().st_size < 1000:
            raise RuntimeError(f"PyMOL output is missing or unexpectedly small: {out_path}")
        print(f"[OK] {tag}: chain {chain}, {expected}{resi}{target}; saved {out_path}")

    print("[DONE] All single-point mutants were built independently from WT.")


main()
