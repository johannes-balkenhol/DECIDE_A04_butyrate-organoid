"""a04_common.py -- shared loaders/helpers for the DECIDE A04 downstream notebooks.
Import in every notebook:  import a04_common as ac
Keeps the (tested) featureCounts parsing + design in ONE place so notebooks stay thin.
"""
import os, re
from pathlib import Path
import numpy as np
import pandas as pd

# ---- paths ------------------------------------------------------------------
def project_root():
    if "A04_PROJ" in os.environ:
        return Path(os.environ["A04_PROJ"]).resolve()
    # notebooks live in notebooks/ -> project root is parent
    here = Path.cwd()
    for p in [here, *here.parents]:
        if (p/"configs/samples.tsv").exists() or (p/"results/tables/counts_raw.tsv").exists():
            return p
    return here

PROJ    = project_root()
COUNTS  = PROJ/"results/tables/counts_raw.tsv"
SAMPLES = PROJ/"configs/samples.tsv"
GTF     = PROJ/"data/refs/human_gencode/gencode.v46.primary_assembly.annotation.gtf"
BAMDIR  = PROJ/"data/processed/bam"
TAB     = PROJ/"results/tables";  FIG = PROJ/"results/figures"
TAB.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)

GROUP_ORDER = ["Mock_0","Mock_5","Mock_10","Cd_0","Cd_5","Cd_10"]
GROUP_COLORS = {"Mock_0":"#9ecae1","Mock_5":"#4292c6","Mock_10":"#08519c",
                "Cd_0":"#fcae91","Cd_5":"#fb6a4a","Cd_10":"#cb181d"}

CONTRASTS = {
    "Cd_vs_Mock_0mM":    ["group","Cd_0","Mock_0"],
    "Cd_vs_Mock_5mM":    ["group","Cd_5","Mock_5"],
    "Cd_vs_Mock_10mM":   ["group","Cd_10","Mock_10"],
    "Butyrate_10v0_Mock":["group","Mock_10","Mock_0"],
    "Butyrate_5v0_Mock": ["group","Mock_5","Mock_0"],
    "Butyrate_10v0_Cd":  ["group","Cd_10","Cd_0"],
    "Butyrate_5v0_Cd":   ["group","Cd_5","Cd_0"],
}
PANELS = {
    "ER chaperones":     ["HSPA5","CALR","WFS1","HYOU1","MANF","PDIA4"],
    "UPR sensors":       ["XBP1","ATF6","ATF4","EIF2AK3"],
    "Regulators":        ["AR","TMPRSS2","BHLHA15"],
    "NF-kB / chemokine": ["NFKB1","RELA","RELB","IRAK2","IKBKE","CXCL1","CXCL2","CXCL3","CXCL8"],
}

# ---- loaders ----------------------------------------------------------------
def load_counts():
    fc = pd.read_csv(COUNTS, sep="\t", comment="#")
    meta_cols = ["Geneid","Chr","Start","End","Strand","Length"]
    counts = fc.drop(columns=[c for c in meta_cols if c in fc.columns]).copy()
    counts.index = fc["Geneid"].values
    counts.columns = [re.sub(r"_Aligned\.sortedByCoord\.out\.bam$","",os.path.basename(c))
                      for c in counts.columns]
    counts.index = [g.split(".")[0] for g in counts.index]
    counts = counts.groupby(level=0).sum()
    return counts

def load_meta(order_like=None):
    sids = list(pd.read_csv(SAMPLES, sep="\t", comment="#")["sample_id"]) \
           if SAMPLES.exists() else (list(order_like) if order_like is not None else [])
    rx = re.compile(r"^(M|Cd)(\d+)R(\d+)$"); rows=[]
    for s in sids:
        m = rx.match(s)
        if not m: continue
        inf = "Mock" if m.group(1)=="M" else "Cd"
        rows.append(dict(sample_id=s, infection=inf, butyrate=int(m.group(2)),
                         replicate=int(m.group(3)), group=f"{inf}_{m.group(2)}"))
    meta = pd.DataFrame(rows).set_index("sample_id")
    meta["group"] = pd.Categorical(meta["group"], categories=GROUP_ORDER, ordered=True)
    meta["infection"] = pd.Categorical(meta["infection"], categories=["Mock","Cd"], ordered=True)
    return meta

def load_symbol_map():
    id2sym = {}
    if not GTF.exists(): return id2sym, {}
    gid = re.compile(r'gene_id "([^"]+)"'); gnm = re.compile(r'gene_name "([^"]+)"')
    with open(GTF) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.split("\t")
            if len(f) > 2 and f[2] == "gene":
                a, b = gid.search(ln), gnm.search(ln)
                if a and b: id2sym[a.group(1).split(".")[0]] = b.group(1)
    return id2sym, {v:k for k,v in id2sym.items()}

def prefilter(counts, min_count=10, min_samples=3):
    return counts.loc[(counts >= min_count).sum(axis=1) >= min_samples]

def lognorm_cpm(counts_f):
    cpm = counts_f / counts_f.sum(0) * 1e6
    return np.log1p(cpm)            # genes x samples

# ---- pydeseq2 (API-version tolerant) ----------------------------------------
def build_dds(counts_f, meta):
    from pydeseq2.dds import DeseqDataSet
    try:
        dds = DeseqDataSet(counts=counts_f.T, metadata=meta, design_factors="group", quiet=True)
    except TypeError:
        dds = DeseqDataSet(counts=counts_f.T, metadata=meta, design="~group", quiet=True)
    dds.deseq2()
    return dds

def run_contrast(dds, contrast, id2sym=None):
    from pydeseq2.ds import DeseqStats
    st = DeseqStats(dds, contrast=contrast, quiet=True); st.summary()
    r = st.results_df.copy()
    if id2sym is not None: r["symbol"] = [id2sym.get(g, g) for g in r.index]
    return r.sort_values("padj")

def parse_star_logs():
    rows = []
    for f in sorted(BAMDIR.glob("*_Log.final.out")):
        sid = f.name.replace("_Log.final.out","")
        d = {"sample_id": sid}
        for ln in open(f):
            if "|" in ln:
                k, v = ln.split("|"); d[k.strip()] = v.strip()
        rows.append(d)
    return pd.DataFrame(rows).set_index("sample_id") if rows else pd.DataFrame()
