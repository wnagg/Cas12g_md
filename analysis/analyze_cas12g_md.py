#!/usr/bin/env python3
"""Compact analysis for Cas12g WT/mutant GROMACS trajectories."""
from __future__ import annotations
import argparse, csv, subprocess, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MUTATIONS = {"F616N": (616,"PHE","ASN"), "R184S": (184,"ARG","SER"),
             "V623Q": (623,"VAL","GLN"), "V623S": (623,"VAL","SER")}
PAIRS = [
 ("R661_R224", "com of (resnr 661 and (name NH1 or name NH2)) plus resnr 224 and name CZ"),
 ("D513_E655", "com of (resnr 513 and (name OD1 or name OD2)) plus com of (resnr 655 and (name OE1 or name OE2))"),
 ("D513_K728", "com of (resnr 513 and (name OD1 or name OD2)) plus resnr 728 and name NZ"),
]
COLORS={"WT":"#0072B2","F616N":"#D55E00","R184S":"#009E73","V623Q":"#CC79A7","V623S":"#E69F00"}

def read_xvg(path):
    rows=[]
    for line in Path(path).read_text(errors="replace").splitlines():
        line=line.strip()
        if not line or line.startswith(("#","@")): continue
        try: rows.append([float(v) for v in line.split()])
        except ValueError: pass
    if not rows: raise RuntimeError(f"No numeric data in {path}")
    return np.asarray(rows,float)

def gmx_run(gmx,args,log,input_text=None):
    cmd=[gmx,*map(str,args)]
    with open(log,"a") as f:
        f.write("\n$ "+" ".join(cmd)+"\n"); f.flush()
        p=subprocess.run(cmd,input=input_text,text=True,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode: raise RuntimeError(f"GROMACS failed ({p.returncode}): {' '.join(cmd)}; see {log}")

def pdb_read(path,chain):
    residues={}; atoms={}; ca={}
    for l in Path(path).read_text(errors="replace").splitlines():
        if not l.startswith("ATOM  ") or len(l)<54 or l[16] not in (" ","A"): continue
        ch=l[21].strip() or "_"
        if ch!=(chain.strip() or "_"): continue
        resi=l[22:26].strip(); ins=l[26].strip(); resn=l[17:20].strip().upper(); atom=l[12:16].strip()
        try: xyz=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        except ValueError: continue
        key=(ch,resi,ins); residues[key]=resn; atoms[(*key,atom)]=xyz
        if atom=="CA": ca[key]=xyz
    if not ca: raise RuntimeError(f"No C-alpha atoms for chain {chain} in {path}")
    return residues,atoms,ca

def fit_rmsd(ref,mob):
    a=ref-ref.mean(0); b=mob-mob.mean(0); u,_,vt=np.linalg.svd(b.T@a); rot=u@vt
    if np.linalg.det(rot)<0: u[:,-1]*=-1; rot=u@vt
    fitted=(mob-mob.mean(0))@rot+ref.mean(0)
    return float(np.sqrt(np.mean(np.sum((fitted-ref)**2,axis=1))))

def initial_qc(structures,outdir,variants,chain):
    outdir.mkdir(parents=True,exist_ok=True)
    wr,wa,wca=pdb_read(structures/"WT.pdb",chain); rows=[]
    for tag in variants:
        p=structures/f"{tag}.pdb"
        if not p.is_file(): raise FileNotFoundError(p)
        r,a,ca=pdb_read(p,chain); pos=expected=target=""
        if tag in MUTATIONS:
            pos,expected,target=MUTATIONS[tag]
            key=(chain,str(pos),"")
            if wr.get(key)!=expected or r.get(key)!=target:
                raise RuntimeError(f"Mutation identity QC failed for {tag}: WT={wr.get(key)}, mutant={r.get(key)}")
        common=sorted(set(wca)&set(ca))
        if len(common)<10: raise RuntimeError(f"Too few matched C-alpha atoms for {tag}")
        rms=fit_rmsd(np.array([wca[k] for k in common]),np.array([ca[k] for k in common]))
        bkeys=sorted(k for k in set(wa)&set(a) if k[3] in ("N","CA","C","O","OXT"))
        br=np.array([wa[k] for k in bkeys]); bm=np.array([a[k] for k in bkeys])
        # Use C-alpha fit to remove rigid-body differences, then quantify shared backbone.
        ca_ref=np.array([wca[k] for k in common]); ca_mob=np.array([ca[k] for k in common])
        c1=ca_ref.mean(0); c2=ca_mob.mean(0); u,_,vt=np.linalg.svd((ca_mob-c2).T@(ca_ref-c1)); rot=u@vt
        if np.linalg.det(rot)<0: u[:,-1]*=-1; rot=u@vt
        bfit=(bm-c2)@rot+c1
        brms=float(np.sqrt(np.mean(np.sum((bfit-br)**2,axis=1))))
        rows.append({"variant":tag,"mutation":f"{expected}{pos}{target}" if pos else "WT","chain":chain,
                     "CA_atoms_compared":len(common),"backbone_atoms_compared":len(bkeys),
                     "global_CA_RMSD_A_after_fit":rms,"backbone_RMSD_A_after_CA_fit":brms,
                     "interpretation":"Construction QC; near-zero backbone RMSD is expected for side-chain mutagenesis"})
    pd.DataFrame(rows).to_csv(outdir/"initial_models_qc.csv",index=False)
    print(f"[OK] {outdir/'initial_models_qc.csv'}")

def analyze(workdir,tag,gmx,deffnm,begin,temp,bins):
    workdir=Path(workdir).resolve(); out=workdir/"analysis"; out.mkdir(parents=True,exist_ok=True)
    tpr=workdir/f"{deffnm}.tpr"; xtc=workdir/f"{deffnm}.xtc"
    if not tpr.is_file() or not xtc.is_file(): raise FileNotFoundError(f"Need production files {tpr} and {xtc}")
    log=out/"gromacs_analysis.log"; log.write_text(f"System={tag}\nTPR={tpr}\nXTC={xtc}\n")
    with tempfile.TemporaryDirectory(prefix=".tmp_",dir=out) as td:
        td=Path(td); pbc=td/"pbc.xtc"; ref=td/"initial.gro"; fit=td/"fit.xtc"
        rms=td/"rmsd.xvg"; rg=td/"rg.xvg"; ds=td/"dist.xvg"; rf=td/"rmsf.xvg"
        gmx_run(gmx,["trjconv","-s",tpr,"-f",xtc,"-o",pbc,"-pbc","mol","-center"],log,"Protein\nSystem\n")
        gmx_run(gmx,["trjconv","-s",tpr,"-f",pbc,"-o",ref,"-dump","0"],log,"System\n")
        gmx_run(gmx,["trjconv","-s",ref,"-f",pbc,"-o",fit,"-fit","rot+trans"],log,"Backbone\nSystem\n")
        gmx_run(gmx,["rms","-s",ref,"-f",fit,"-o",rms,"-tu","ns"],log,"Backbone\nBackbone\n")
        gmx_run(gmx,["gyrate","-s",tpr,"-f",fit,"-o",rg,"-tu","ns"],log,"Protein\n")
        gmx_run(gmx,["distance","-s",tpr,"-f",fit,"-select","; ".join(x[1] for x in PAIRS),"-oall",ds,"-tu","ns"],log)
        r=read_xvg(rms); g=read_xvg(rg); d=read_xvg(ds)
        if d.shape[1]<4: raise RuntimeError(f"Expected 3 distance columns, got {d.shape[1]}; see {log}")
        time=r[:,0]; duration=float(time.max()); start=min(begin,max(0,duration*0.2),max(0,duration-1e-6))
        gmx_run(gmx,["rmsf","-s",ref,"-f",fit,"-o",rf,"-res","-b",f"{start*1000:.6f}"],log,"C-alpha\n")
        f=read_xvg(rf); rgval=np.interp(time,g[:,0],g[:,1]); dist=[np.interp(time,d[:,0],d[:,i]) for i in (1,2,3)]
        metrics=pd.DataFrame({"time_ns":time,"RMSD_nm":r[:,1],"Rg_nm":rgval,
            "R661_R224_nm":dist[0],"D513_E655_nm":dist[1],"D513_K728_nm":dist[2]})
        metrics.to_csv(out/f"metrics_{tag}.csv",index=False)
        rmsf=pd.DataFrame({"residue":f[:,0],"RMSF_nm":f[:,1]}); rmsf.to_csv(out/f"rmsf_{tag}.csv",index=False)
        sel=metrics[metrics.time_ns>=start]
        nb=max(5,min(int(bins),int(np.sqrt(max(1,len(sel))/2))))
        H,xe,ye=np.histogram2d(sel.RMSD_nm,sel.Rg_nm,bins=nb); prob=H/H.sum()
        fe=np.full_like(H,np.nan); occ=H>0; fe[occ]=-8.314462618e-3*temp*np.log(prob[occ]); fe[occ]-=np.nanmin(fe[occ])
        pd.DataFrame([{"RMSD_nm":(xe[i]+xe[i+1])/2,"Rg_nm":(ye[j]+ye[j+1])/2,"count":int(H[i,j]),
            "free_energy_kJmol":fe[i,j] if occ[i,j] else np.nan} for i in range(nb) for j in range(nb)]).to_csv(out/f"FEL_RMSD_Rg_{tag}.csv",index=False)
        fig,ax=plt.subplots(2,3,figsize=(13,8),constrained_layout=True); ax=ax.ravel(); color=COLORS.get(tag)
        for a,col,ylabel,title in [(ax[0],"RMSD_nm","RMSD (nm)","Backbone RMSD"),(ax[1],"Rg_nm","Rg (nm)","Radius of gyration")]:
            a.plot(time,metrics[col],lw=.8,color=color); a.set(xlabel="Time (ns)",ylabel=ylabel,title=title)
        ax[2].plot(rmsf.residue,rmsf.RMSF_nm,lw=.8,color=color); ax[2].set(xlabel="Residue",ylabel="Cα RMSF (nm)",title=f"RMSF after {start:g} ns")
        mesh=ax[3].pcolormesh(xe,ye,np.ma.masked_invalid(fe.T),shading="auto",cmap="viridis")
        fig.colorbar(mesh,ax=ax[3],label="Relative free energy (kJ/mol)"); ax[3].set(xlabel="RMSD (nm)",ylabel="Rg (nm)",title="RMSD–Rg FEL")
        for a,col,title in zip(ax[4:],("R661_R224_nm","D513_E655_nm","D513_K728_nm"),("R661–R224","D513–E655","D513–K728")):
            a.plot(time,metrics[col],lw=.8,color=color); a.set(xlabel="Time (ns)",ylabel="Distance (nm)",title=title)
        fig.suptitle(f"Cas12g {tag}",fontsize=14); fig.savefig(out/f"summary_{tag}.pdf"); fig.savefig(out/f"summary_{tag}.png",dpi=300); plt.close(fig)
        (out/f"summary_{tag}.txt").write_text(f"System: {tag}\nDuration_ns: {duration:.3f}\nFrames: {len(metrics)}\nRMSD reference: first production frame\nAnalysis start ns: {start:.3f}\nFEL bins: {nb}\n")
    print(f"[OK] {tag} analysis complete: {out}")

def compare(root,out,variants,begin,temp):
    out.mkdir(parents=True,exist_ok=True); data={}
    for tag in variants:
        p=Path(root)/"runs"/tag/"analysis"/f"metrics_{tag}.csv"
        if not p.is_file(): raise FileNotFoundError(f"Missing {p}; run MD and analysis first")
        data[tag]=pd.read_csv(p)
    fig,axs=plt.subplots(2,3,figsize=(14,8),constrained_layout=True); axs=axs.ravel()
    specs=[("RMSD_nm","RMSD (nm)"),("Rg_nm","Rg (nm)"),("R661_R224_nm","R661–R224 (nm)"),("D513_E655_nm","D513–E655 (nm)"),("D513_K728_nm","D513–K728 (nm)")]
    for ax,(col,label) in zip(axs,specs):
        for tag,df in data.items(): ax.plot(df.time_ns,df[col],lw=.7,label=tag,color=COLORS.get(tag))
        ax.set(xlabel="Time (ns)",ylabel=label,title=label)
    for tag in variants:
        p=Path(root)/"runs"/tag/"analysis"/f"rmsf_{tag}.csv"; df=pd.read_csv(p)
        axs[5].plot(df.residue,df.RMSF_nm,lw=.8,label=tag,color=COLORS.get(tag))
    axs[5].set(xlabel="Residue",ylabel="Cα RMSF (nm)",title="Cα RMSF")
    handles,labels=axs[0].get_legend_handles_labels(); fig.legend(handles,labels,loc="upper center",ncol=len(variants),frameon=False)
    fig.savefig(out/"comparison_timeseries.pdf",bbox_inches="tight"); fig.savefig(out/"comparison_timeseries.png",dpi=300,bbox_inches="tight"); plt.close(fig)
    # Common axes and color scale for a fair visual comparison; empty bins remain masked.
    xs=np.concatenate([df.RMSD_nm.to_numpy() for df in data.values()]); ys=np.concatenate([df.Rg_nm.to_numpy() for df in data.values()])
    xr=(float(xs.min()),float(xs.max())); yr=(float(ys.min()),float(ys.max())); nb=max(5,min(30,int(np.sqrt(min(len(df) for df in data.values())/2))))
    grids=[]
    for tag,df in data.items():
        duration=float(df.time_ns.max()); start=min(float(begin),max(0.0,duration*0.2),max(0.0,duration-1e-6))
        df=df[df.time_ns>=start]
        if len(df)<10: raise RuntimeError(f'Too few frames after analysis start for {tag}')
        H,xe,ye=np.histogram2d(df.RMSD_nm,df.Rg_nm,bins=nb,range=[xr,yr]); P=H/H.sum(); fe=np.full_like(H,np.nan); occ=H>0
        fe[occ]=-8.314462618e-3*float(temp)*np.log(P[occ]); fe[occ]-=np.nanmin(fe[occ]); grids.append((tag,xe,ye,fe))
    vals=np.concatenate([g[3][np.isfinite(g[3])] for g in grids]); vmax=max(5,float(np.percentile(vals,95)))
    fig,axs=plt.subplots(2,3,figsize=(13,8),constrained_layout=True); axs=axs.ravel(); mesh=None
    for ax,(tag,xe,ye,fe) in zip(axs,grids):
        mesh=ax.pcolormesh(xe,ye,np.ma.masked_invalid(fe.T),shading="auto",cmap="viridis",vmin=0,vmax=vmax); ax.set(title=tag,xlabel="RMSD (nm)",ylabel="Rg (nm)")
    for ax in axs[len(grids):]: ax.axis("off")
    fig.colorbar(mesh,ax=axs[:len(grids)],label="Relative free energy (kJ/mol)",shrink=.85)
    fig.savefig(out/"comparison_FEL.pdf",bbox_inches="tight"); fig.savefig(out/"comparison_FEL.png",dpi=300,bbox_inches="tight"); plt.close(fig)
    print(f"[OK] Comparison figures: {out}")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=["initial","run","compare"],required=True); p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--structures-dir",type=Path); p.add_argument("--variants",nargs="+",default=["WT","F616N","R184S","V623Q","V623S"]); p.add_argument("--chain",default="A")
    p.add_argument("--workdir",type=Path); p.add_argument("--root",type=Path); p.add_argument("--tag"); p.add_argument("--gmx",default="gmx"); p.add_argument("--deffnm",default="md500")
    p.add_argument("--begin-ns",type=float,default=20); p.add_argument("--temp-k",type=float,default=300); p.add_argument("--fel-bins",type=int,default=30)
    a=p.parse_args()
    if a.mode=="initial":
        if not a.structures_dir: p.error("--structures-dir required")
        initial_qc(a.structures_dir,a.outdir,a.variants,a.chain)
    elif a.mode=="run":
        if not a.workdir or not a.tag: p.error("--workdir and --tag required")
        analyze(a.workdir,a.tag,a.gmx,a.deffnm,a.begin_ns,a.temp_k,a.fel_bins)
    else:
        if not a.root: p.error("--root required")
        compare(a.root,a.outdir,a.variants,a.begin_ns,a.temp_k)
if __name__=="__main__": main()
