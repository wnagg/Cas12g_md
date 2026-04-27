#!/usr/bin/env bash
# Reproducible GROMACS workflow for protein-only Cas12g MD.
# Usage: bash scripts/run_md_new_protein.sh Cas12g_predict_F616N.env

set -Eeuo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash $0 /path/to/config.env"
  exit 1
fi

CONFIG=$1
source "$CONFIG"

mkdir -p "$WORKDIR" "$WORKDIR/logs" "$WORKDIR/.done"
cd "$WORKDIR"

stage(){ echo "$1" | tee CURRENT_STAGE; }
done_stage(){ touch ".done/$1.done"; }
has_done(){ [[ -f ".done/$1.done" ]]; }

ps_to_steps(){ python3 - <<PY
ps=float('$1'); dt=float('${DT_PS:-0.002}'); print(int(round(ps/dt)))
PY
}
ns_to_steps(){ python3 - <<PY
ns=float('$1'); dt=float('${DT_PS:-0.002}'); print(int(round(ns*1000/dt)))
PY
}
interval_to_steps(){ python3 - <<PY
ps=float('$1'); dt=float('${DT_PS:-0.002}'); print(max(1,int(round(ps/dt))))
PY
}

cat > md_monitor.sh <<'EOF'
#!/usr/bin/env bash
echo "Time: $(date)"
echo "Current stage: $(cat CURRENT_STAGE 2>/dev/null || echo unknown)"
echo
ls .done 2>/dev/null | sed 's/.done$//' || true
echo
echo "Recent file sizes:"
ls -lh *.gro *.xtc *.edr *.log *.tpr *.cpt 2>/dev/null | awk '{print $5, $9}' | tail -n 25
EOF
chmod +x md_monitor.sh

if ! has_done 00_write_mdp; then
  stage 00_write_mdp
  cat > ions.mdp <<EOF
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 5000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF
  cat > minim.mdp <<EOF
integrator  = steep
nsteps      = 50000
emtol       = 1000.0
emstep      = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
constraints = none
pbc         = xyz
EOF
  cat > nvt.mdp <<EOF
define       = -DPOSRES_CA
integrator   = md
dt           = ${DT_PS:-0.002}
nsteps       = $(ps_to_steps ${NVT_PS:-100})
continuation = no
gen_vel      = yes
gen_temp     = ${TEMP_K:-300}
gen_seed     = -1
cutoff-scheme = Verlet
nstlist      = 20
rlist        = 1.0
coulombtype  = PME
rcoulomb     = 1.0
vdwtype      = Cut-off
rvdw         = 1.0
pbc          = xyz
constraints  = h-bonds
constraint_algorithm = lincs
lincs_iter   = 1
lincs_order  = 4
tcoupl       = V-rescale
tc-grps      = Protein Water_and_ions
tau_t        = 0.1   0.1
ref_t        = ${TEMP_K:-300} ${TEMP_K:-300}
pcoupl       = no
nstxout      = 0
nstvout      = 0
nstenergy    = $(interval_to_steps ${ENERGY_INTERVAL_PS:-10})
nstlog       = $(interval_to_steps ${LOG_INTERVAL_PS:-10})
nstxout-compressed = 500
EOF
  cat > npt.mdp <<EOF
define       = -DPOSRES_CA
integrator   = md
dt           = ${DT_PS:-0.002}
nsteps       = $(ps_to_steps ${NPT_CA_PS:-1000})
continuation = yes
gen_vel      = no
cutoff-scheme = Verlet
nstlist      = 20
rlist        = 1.0
coulombtype  = PME
rcoulomb     = 1.0
vdwtype      = Cut-off
rvdw         = 1.0
pbc          = xyz
constraints  = h-bonds
constraint_algorithm = lincs
lincs_iter   = 1
lincs_order  = 4
tcoupl       = V-rescale
tc-grps      = Protein Water_and_ions
tau_t        = 0.1   0.1
ref_t        = ${TEMP_K:-300} ${TEMP_K:-300}
pcoupl       = C-rescale
pcoupltype   = isotropic
tau_p        = 5.0
ref_p        = ${PRESSURE_BAR:-1.0}
compressibility = 4.5e-5
refcoord_scaling = com
nstxout      = 0
nstvout      = 0
nstenergy    = 1000
nstlog       = 1000
nstxout-compressed = 500
EOF
  cat > npt_500ps_noca.mdp <<EOF
integrator   = md
dt           = ${DT_PS:-0.002}
nsteps       = $(ps_to_steps ${NPT_NOCA_PS:-500})
continuation = yes
gen_vel      = no
cutoff-scheme = Verlet
nstlist      = 20
rlist        = 1.0
coulombtype  = PME
rcoulomb     = 1.0
vdwtype      = Cut-off
rvdw         = 1.0
pbc          = xyz
constraints  = h-bonds
constraint_algorithm = lincs
lincs_iter   = 1
lincs_order  = 4
tcoupl       = V-rescale
tc-grps      = Protein Water_and_ions
tau_t        = 0.1   0.1
ref_t        = ${TEMP_K:-300} ${TEMP_K:-300}
pcoupl       = C-rescale
pcoupltype   = isotropic
tau_p        = 5.0
ref_p        = ${PRESSURE_BAR:-1.0}
compressibility = 4.5e-5
nstxout      = 0
nstvout      = 0
nstenergy    = 1000
nstlog       = 1000
nstxout-compressed = 500
EOF
  cat > md_500ns.mdp <<EOF
integrator   = md
dt           = ${DT_PS:-0.002}
nsteps       = $(ns_to_steps ${MD_NS:-500})
continuation = yes
gen_vel      = no
cutoff-scheme = Verlet
nstlist      = 100
rlist        = 1.0
coulombtype  = PME
rcoulomb     = 1.0
vdwtype      = Cut-off
rvdw         = 1.0
DispCorr     = EnerPres
pbc          = xyz
constraints  = h-bonds
constraint_algorithm = lincs
lincs_iter   = 1
lincs_order  = 4
tcoupl       = V-rescale
tc-grps      = Protein Water_and_ions
tau_t        = 0.1   0.1
ref_t        = ${TEMP_K:-300} ${TEMP_K:-300}
pcoupl       = C-rescale
pcoupltype   = isotropic
tau_p        = 5.0
ref_p        = ${PRESSURE_BAR:-1.0}
compressibility = 4.5e-5
nstxout      = 0
nstvout      = 0
nstfout      = 0
nstenergy    = $(interval_to_steps ${ENERGY_INTERVAL_PS:-10})
nstlog       = $(interval_to_steps ${LOG_INTERVAL_PS:-10})
nstxout-compressed = $(interval_to_steps ${XTC_INTERVAL_PS:-100})
EOF
  done_stage 00_write_mdp
fi

if ! has_done 01_prepare_input; then
  stage 01_prepare_input
  if [[ "${PREPARE_PROTEIN_ONLY:-yes}" == "yes" ]]; then
    grep -E '^(ATOM  |TER   |END   )' "$INPUT_PDB" > input_protein.pdb
  else
    cp "$INPUT_PDB" input_protein.pdb
  fi
  done_stage 01_prepare_input
fi

if ! has_done 02_topology_solvate_ions; then
  stage 02_topology_solvate_ions
  "$GMX" pdb2gmx -f input_protein.pdb -o processed.gro -p topol.top -ff ${FORCEFIELD:-amber99sb-ildn} -water ${WATER:-tip3p} -ignh 2>&1 | tee logs/pdb2gmx.out
  "$GMX" editconf -f processed.gro -o newbox.gro -c -d ${BOX_DIST_NM:-1.0} -bt ${BOX_TYPE:-dodecahedron} 2>&1 | tee logs/editconf.out
  "$GMX" solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tee logs/solvate.out
  "$GMX" grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_ions.out
  printf "SOL\n" | "$GMX" genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc ${SALT_MOLAR:-0.15} 2>&1 | tee logs/genion.out
  done_stage 02_topology_solvate_ions
fi

if ! has_done 03_em; then
  stage 03_em
  "$GMX" grompp -f minim.mdp -c solvated_ions.gro -p topol.top -o em.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_em.out
  "$GMX" mdrun -deffnm em -v -nb cpu -pme cpu 2>&1 | tee logs/em.out
  done_stage 03_em
fi

if ! has_done 04_make_ca_restraint; then
  stage 04_make_ca_restraint
  printf "3\nq\n" | "$GMX" make_ndx -f em.gro -o index_ca.ndx 2>&1 | tee logs/make_ndx_ca.out
  printf "3\n" | "$GMX" genrestr -f em.gro -n index_ca.ndx -o posre_Ca.itp -fc ${CA_FC:-1000} ${CA_FC:-1000} ${CA_FC:-1000} 2>&1 | tee logs/genrestr_ca.out
  python3 - <<'PY'
from pathlib import Path
p=Path('topol.top')
s=p.read_text()
inc='#ifdef POSRES_CA\n#include "posre_Ca.itp"\n#endif\n\n'
anchor='; Include water topology'
if 'posre_Ca.itp' not in s:
    s=s.replace(anchor, inc+anchor)
    p.write_text(s)
PY
  done_stage 04_make_ca_restraint
fi

if ! has_done 05_nvt; then
  stage 05_nvt
  "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_nvt.out
  "$GMX" mdrun -deffnm nvt -v -cpt ${CHECKPOINT_MIN:-15} ${MDRUN_CPU_OPTS:-"-nb cpu -pme cpu"} 2>&1 | tee logs/nvt.out
  done_stage 05_nvt
fi

if ! has_done 06_npt_ca; then
  stage 06_npt_ca
  "$GMX" grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_npt.out
  "$GMX" mdrun -deffnm npt -v -cpt ${CHECKPOINT_MIN:-15} ${MDRUN_CPU_OPTS:-"-nb cpu -pme cpu"} 2>&1 | tee logs/npt.out
  done_stage 06_npt_ca
fi

if ! has_done 07_npt_noca; then
  stage 07_npt_noca
  "$GMX" grompp -f npt_500ps_noca.mdp -c npt.gro -t npt.cpt -p topol.top -o npt_500ps_noca.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_npt_noca.out
  "$GMX" mdrun -deffnm npt_500ps_noca -v -cpt ${CHECKPOINT_MIN:-15} ${MDRUN_CPU_OPTS:-"-nb cpu -pme cpu"} 2>&1 | tee logs/npt_noca.out
  done_stage 07_npt_noca
fi

if ! has_done 09_prepare_production_start; then
  stage 09_prepare_production_start
  printf "1\n0\n" | "$GMX" trjconv -s npt_500ps_noca.tpr -f npt_500ps_noca.gro -o npt_500ps_noca_center.gro -pbc mol -center 2>&1 | tee logs/center_start.out
  done_stage 09_prepare_production_start
fi

if ! has_done 10_md500; then
  stage 10_md500
  "$GMX" grompp -f md_500ns.mdp -c npt_500ps_noca_center.gro -t npt_500ps_noca.cpt -p topol.top -o md500.tpr -maxwarn ${MAXWARN:-1} 2>&1 | tee logs/grompp_md500.out
  "$GMX" mdrun -deffnm md500 -v -cpt ${CHECKPOINT_MIN:-15} ${MDRUN_GPU_OPTS:-"-nb gpu -pme gpu -bonded gpu -update gpu"} 2>&1 | tee logs/md500.out
  done_stage 10_md500
fi

stage ALL_DONE
echo "ALL DONE: $WORKDIR"
