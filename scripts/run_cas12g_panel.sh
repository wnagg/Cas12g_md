#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"; REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE="$REPO_ROOT/scripts/run_md_new_protein.sh"; ANALYZER="$REPO_ROOT/analysis/analyze_cas12g_md.py"
CONFIG="$REPO_ROOT/env/Cas12g_panel.env"; MODE=all
if [[ $# -ge 1 && "$1" != --* ]]; then CONFIG="$1"; shift; fi
if [[ $# -ge 1 ]]; then MODE="$1"; fi
if [[ "$MODE" == -h || "$MODE" == --help ]]; then echo "Usage: $0 [config.env] [all|--build-only|--run-only|--analyze-only]"; exit 0; fi
case "$MODE" in all|--build-only|--run-only|--analyze-only) ;; *) echo "Unknown mode $MODE" >&2; exit 2;; esac
[[ -f "$CONFIG" ]] || { echo "Copy env/Cas12g_panel.env.example to env/Cas12g_panel.env and edit paths." >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"
for var in GMX PYMOL PYTHON WT_PDB CHAIN WORKDIR VARIANTS REBUILD_MUTANTS FORCEFIELD WATER BOX_TYPE BOX_DIST_NM SALT_MOLAR TEMP_K PRESSURE_BAR PREPARE_PROTEIN_ONLY NVT_PS NPT_CA_PS NPT_NOCA_PS MD_NS DT_PS XTC_INTERVAL_PS ENERGY_INTERVAL_PS LOG_INTERVAL_PS CHECKPOINT_MIN MD_DEFFNM CA_FC USE_GATE_CA GATE_RESID GATE_FC MDRUN_CPU_OPTS MDRUN_GPU_OPTS MAXWARN ANALYSIS_BEGIN_NS FEL_BINS CLEAN_INTERMEDIATES; do
  [[ -n "$(eval echo \$$var)" ]] || { echo "Missing config variable: $var" >&2; exit 1; }
done
[[ -f "$WT_PDB" ]] || { echo "WT_PDB not found: $WT_PDB" >&2; exit 1; }
WT_PDB="$(cd "$(dirname "$WT_PDB")" && printf '%s/%s' "$PWD" "$(basename "$WT_PDB")")"
WORKDIR="$(mkdir -p "$WORKDIR" && cd "$WORKDIR" && pwd)"
STRUCT_DIR="$WORKDIR/structures"; RUN_ROOT="$WORKDIR/runs"; ANALYSIS_DIR="$WORKDIR/analysis"; CONFIG_DIR="$WORKDIR/.configs"
mkdir -p "$STRUCT_DIR" "$RUN_ROOT" "$ANALYSIS_DIR" "$CONFIG_DIR"
[[ -x "$GMX" ]] || { echo "GROMACS executable unavailable: $GMX" >&2; exit 1; }
command -v "$PYTHON" >/dev/null || { echo "Python unavailable: $PYTHON" >&2; exit 1; }
"$PYTHON" -c 'import numpy,pandas,matplotlib'
for tag in $VARIANTS; do case "$tag" in WT|F616N|R184S|V623Q|V623S) ;; *) echo "Unsupported variant: $tag" >&2; exit 1;; esac; done
existing_tpr="$(find "$RUN_ROOT" -type f -name "$MD_DEFFNM.tpr" -print -quit 2>/dev/null || true)"
build_needed=no
for tag in WT F616N R184S V623Q V623S; do [[ -s "$STRUCT_DIR/$tag.pdb" ]] || build_needed=yes; done
if [[ "$REBUILD_MUTANTS" == yes ]]; then build_needed=yes; fi
if [[ "$MODE" != --run-only && "$MODE" != --analyze-only && "$build_needed" == yes ]]; then
 [[ -z "$existing_tpr" ]] || { echo "Existing production TPR found; use a new WORKDIR before rebuilding." >&2; exit 1; }
 command -v "$PYMOL" >/dev/null || { echo "PyMOL command unavailable: $PYMOL" >&2; exit 1; }
 export CAS12G_WT_PDB="$WT_PDB" CAS12G_STRUCT_DIR="$STRUCT_DIR" CAS12G_CHAIN="$CHAIN"
 (cd "$REPO_ROOT" && "$PYMOL" -cq -d "run scripts/build_cas12g_mutants.py" -d "quit") > "$WORKDIR/build_mutants.log" 2>&1 || { tail -n 60 "$WORKDIR/build_mutants.log" >&2; exit 1; }
 tail -n 20 "$WORKDIR/build_mutants.log"
fi
if [[ "$MODE" == --build-only ]]; then "$PYTHON" "$ANALYZER" --mode initial --structures-dir "$STRUCT_DIR" --outdir "$ANALYSIS_DIR" --variants $VARIANTS --chain "$CHAIN"; exit 0; fi
write_kv(){ printf '%s=%q\n' "$1" "$2"; }
write_run_config(){ tag="$1"; pdb="$2"; wd="$3"; out="$4"; {
 write_kv GMX "$GMX"; write_kv INPUT_PDB "$pdb"; write_kv WORKDIR "$wd"; write_kv SYSTEM_NAME "Cas12g_$tag"
 write_kv FORCEFIELD "$FORCEFIELD"; write_kv WATER "$WATER"; write_kv BOX_TYPE "$BOX_TYPE"; write_kv BOX_DIST_NM "$BOX_DIST_NM"
 write_kv SALT_MOLAR "$SALT_MOLAR"; write_kv TEMP_K "$TEMP_K"; write_kv PRESSURE_BAR "$PRESSURE_BAR"; write_kv PREPARE_PROTEIN_ONLY "$PREPARE_PROTEIN_ONLY"
 write_kv NVT_PS "$NVT_PS"; write_kv NPT_CA_PS "$NPT_CA_PS"; write_kv NPT_NOCA_PS "$NPT_NOCA_PS"; write_kv MD_NS "$MD_NS"; write_kv DT_PS "$DT_PS"
 write_kv XTC_INTERVAL_PS "$XTC_INTERVAL_PS"; write_kv ENERGY_INTERVAL_PS "$ENERGY_INTERVAL_PS"; write_kv LOG_INTERVAL_PS "$LOG_INTERVAL_PS"; write_kv CHECKPOINT_MIN "$CHECKPOINT_MIN"
 write_kv MD_DEFFNM "$MD_DEFFNM"; write_kv CA_FC "$CA_FC"; write_kv USE_GATE_CA "$USE_GATE_CA"; write_kv GATE_RESID "$GATE_RESID"; write_kv GATE_FC "$GATE_FC"
 write_kv MDRUN_CPU_OPTS "$MDRUN_CPU_OPTS"; write_kv MDRUN_GPU_OPTS "$MDRUN_GPU_OPTS"; write_kv MAXWARN "$MAXWARN"; write_kv RUN_LIGHT_ANALYSIS no
 write_kv ANALYSIS_BEGIN_NS "$ANALYSIS_BEGIN_NS"; write_kv CLEAN_EQUIL_XTC yes
} > "$out"; }
"$PYTHON" "$ANALYZER" --mode initial --structures-dir "$STRUCT_DIR" --outdir "$ANALYSIS_DIR" --variants $VARIANTS --chain "$CHAIN"
if [[ "$MODE" == --analyze-only ]]; then
 for tag in $VARIANTS; do "$PYTHON" "$ANALYZER" --mode run --workdir "$RUN_ROOT/$tag" --tag "$tag" --gmx "$GMX" --deffnm "$MD_DEFFNM" --begin-ns "$ANALYSIS_BEGIN_NS" --temp-k "$TEMP_K" --fel-bins "$FEL_BINS" --outdir "$RUN_ROOT/$tag/analysis"; done
else
 for tag in $VARIANTS; do
  pdb="$STRUCT_DIR/$tag.pdb"; wd="$RUN_ROOT/$tag"; cfg="$CONFIG_DIR/$tag.env"; mkdir -p "$wd"
  [[ -s "$pdb" ]] || { echo "Missing $pdb; build structures first." >&2; exit 1; }
  write_run_config "$tag" "$pdb" "$wd" "$cfg"
  echo "=== $tag: $MD_NS ns, dt=$DT_PS ps ==="
  bash "$PIPELINE" "$cfg" 2>&1 | tee "$wd/pipeline.log"
  "$PYTHON" "$ANALYZER" --mode run --workdir "$wd" --tag "$tag" --gmx "$GMX" --deffnm "$MD_DEFFNM" --begin-ns "$ANALYSIS_BEGIN_NS" --temp-k "$TEMP_K" --fel-bins "$FEL_BINS" --outdir "$wd/analysis"
  if [[ "$CLEAN_INTERMEDIATES" == yes ]]; then
   rm -f "$wd"/newbox.gro "$wd"/solvated.gro "$wd"/solvated_ions.gro "$wd"/ions.tpr "$wd"/em.tpr
   rm -f "$wd"/nvt.tpr "$wd"/npt.tpr "$wd"/npt_500ps_noca.tpr "$wd"/nvt.gro "$wd"/npt.gro "$wd"/npt_500ps_noca.gro "$wd"/npt_500ps_noca_center.gro
   rm -f "$wd"/nvt.xtc "$wd"/npt.xtc "$wd"/npt_500ps_noca.xtc "$wd"/em.edr "$wd"/em.log "$wd"/em.cpt "$wd"/nvt.edr "$wd"/nvt.log "$wd"/nvt.cpt "$wd"/npt.edr "$wd"/npt.log "$wd"/npt.cpt "$wd"/npt_500ps_noca.edr "$wd"/npt_500ps_noca.log "$wd"/npt_500ps_noca.cpt
  fi
 done
fi
"$PYTHON" "$ANALYZER" --mode compare --root "$WORKDIR" --outdir "$ANALYSIS_DIR" --variants $VARIANTS
echo "[DONE] $WORKDIR"
