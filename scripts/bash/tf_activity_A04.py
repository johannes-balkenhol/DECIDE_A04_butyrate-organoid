#!/usr/bin/env python3
# =============================================================================
# tf_activity_A04.py — who controls the coupled UPR/chemokine module?
# TF-activity inference with decoupleR + CollecTRI on A04 contrasts.
# Nominates TFs whose activity FALLS as butyrate rises in Cd (candidate shared
# controllers of both the chaperone and chemokine arms), restricted to genes
# actually expressed in A04 (no topology-only nodes -> the A06 AR trap avoided).
#
# Run in an env with: decoupler, pandas, numpy  (pip install decoupler omnipath)
#   A04_PROJ=~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate python scripts/tf_activity_A04.py
# Uses the DE_*.csv stat (Wald 'stat' or signed -log10p) as the gene-level input.
# =============================================================================
import os, sys, re
from pathlib import Path
import numpy as np, pandas as pd

PROJ = Path(os.environ.get("A04_PROJ", Path.cwd()))
TAB  = PROJ/"results/tables"; FIG = PROJ/"results/figures"
OUT  = PROJ/"results/tf_activity"; OUT.mkdir(parents=True, exist_ok=True)
GTF  = PROJ/"data/refs/human_gencode/gencode.v46.primary_assembly.annotation.gtf"

# ---- map Ensembl -> symbol (decoupler/CollecTRI work on symbols) ------------
def gtf_symbol_map(gtf):
    d={}; gid=re.compile(r'gene_id "([^".]+)'); gnm=re.compile(r'gene_name "([^"]+)"')
    if not Path(gtf).exists(): return d
    with open(gtf) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f=ln.split("\t")
            if len(f)>2 and f[2]=="gene":
                a,b=gid.search(ln),gnm.search(ln)
                if a and b: d[a.group(1)]=b.group(1)
    return d
sym = gtf_symbol_map(GTF)

try:
    import decoupler as dc
except ImportError:
    sys.exit("decoupler not installed -> pip install decoupler omnipath")

# ---- CollecTRI regulons (expressed-gene restricted later) -------------------
print(">> fetching CollecTRI...")
net = dc.get_collectri(organism="human", split_complexes=False)

contrasts = ["Butyrate_10v0_Cd","Butyrate_10v0_Mock",
             "Cd_vs_Mock_0mM","Cd_vs_Mock_10mM"]
def load_stat(tag):
    f = TAB/f"DE_{tag}.csv"
    if not f.exists(): print(f"   skip (missing) {f}"); return None
    df = pd.read_csv(f, index_col=0)
    # gene-level signed statistic: prefer Wald 'stat', else signed -log10(padj)
    if "stat" in df.columns: s = df["stat"]
    else:
        s = np.sign(df["log2FoldChange"]) * -np.log10(df["padj"].clip(lower=1e-300))
    s.index = [sym.get(re.sub(r"\..*","",g), g) for g in s.index]
    s = s[~s.index.duplicated()].dropna()
    return s

mats = {t: load_stat(t) for t in contrasts}
mats = {t:s for t,s in mats.items() if s is not None}
mat = pd.DataFrame(mats).T.fillna(0.0)          # contrasts x genes
print(f">> stat matrix: {mat.shape[0]} contrasts x {mat.shape[1]} genes")

# restrict regulons to genes present (expressed) in A04 -> no phantom TFs
net_expr = net[net["target"].isin(mat.columns)]
print(f">> CollecTRI edges kept (expressed targets): {len(net_expr)} / {len(net)}")

# ---- ULM activity -----------------------------------------------------------
acts, pvals = dc.run_ulm(mat=mat, net=net_expr, source="source", target="target",
                         weight="weight", min_n=5, verbose=True)
acts.to_csv(OUT/"tf_activity_scores.csv")
pvals.to_csv(OUT/"tf_activity_pvals.csv")

# ---- the question: which TFs LOSE activity as butyrate rises in Cd? ---------
key = "Butyrate_10v0_Cd"
if key in acts.index:
    a = acts.loc[key].sort_values()
    p = pvals.loc[key]
    res = pd.DataFrame({"activity": a, "pval": p.reindex(a.index)})
    res["fdr"] = (res["pval"]*len(res)).clip(upper=1)   # rough BH-ish
    res.to_csv(OUT/"TF_butyrate_in_Cd_ranked.csv")
    print("\n#### TFs DOWN with butyrate in Cd (candidate module controllers) ####")
    print(res.head(15).round(3).to_string())
    print("\n#### TFs UP with butyrate in Cd ####")
    print(res.tail(8).round(3).to_string())
    # sanity: do the usual inflammatory TFs (NFKB1/RELA) drop? does any chaperone TF drop?
    for tf in ["NFKB1","RELA","RELB","SP1","SREBF2","XBP1","ATF4","ATF6","CREB3L2","NFE2L2","NFYA","CEBPB"]:
        if tf in acts.columns:
            print(f"   {tf:8s} activity@but10v0_Cd = {acts.loc[key,tf]:+.2f}  (p={pvals.loc[key,tf]:.2g})")

# ---- cross-contrast heatmap of top movers ----------------------------------
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    top = acts.loc[key].abs().sort_values(ascending=False).head(25).index if key in acts.index else acts.columns[:25]
    H = acts[top].T
    fig,ax=plt.subplots(figsize=(7, .3*len(top)+1))
    im=ax.imshow(H.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=7)
    fig.colorbar(im, label="TF activity (ULM)"); ax.set_title("TF activity across A04 contrasts")
    plt.tight_layout(); fig.savefig(FIG/"tf_activity_heatmap.png", dpi=130)
    print(f"\n>> wrote {FIG/'tf_activity_heatmap.png'} and {OUT}/TF_butyrate_in_Cd_ranked.csv")
except Exception as e:
    print("plot skipped:", e)
print("\nDONE. The TFs at the TOP of the 'down with butyrate in Cd' list that ALSO")
print("rise with infection (Cd_vs_Mock_0mM) are the candidate shared controllers.")
