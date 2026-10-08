#!/usr/bin/env python3
"""
Generate a FAIR raw-data metadata YAML for a DECIDE project.
Scans data/raw_data/ for fastq files, infers sample IDs, writes
data/meta_data/<project>_raw_metadata.yaml with:
  - dataset-level biological + technical metadata (standardized vocab where possible)
  - per-sample records
  - per-fastq file pointers (relative path) linked to their sample

Usage:
  python make_metadata.py                # auto-detect project from cwd
  python make_metadata.py /path/to/PROJECT
Edit the CONFIG block near the top to match the study, then re-run.
"""
import os, re, sys, hashlib, glob
from pathlib import Path
from datetime import date

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
NAME = ROOT.name
RAW  = ROOT / "data" / "raw_data"
META = ROOT / "data" / "meta_data"; META.mkdir(parents=True, exist_ok=True)

# ---------------- CONFIG: per-project study metadata ----------------
# Filled from what we know; <TODO> marks fields to confirm before publishing.
STUDIES = {
  "DECIDE_A04_CdiffOrganoid_Butyrate": dict(
    title="C. difficile infection of colonic organoids under a butyrate dose-response",
    decide_project="A04",
    organism="Homo sapiens", ncbi_taxon="NCBI:txid9606",
    tissue="colon", tissue_ontology="UBERON:0001155",
    model_system="patient-derived colonic organoid",
    cell_line=None, donor="single donor (<TODO donor ID>)",
    pathogen="Clostridioides difficile", pathogen_taxon="NCBI:txid1496",
    treatment="sodium butyrate", treatment_doses_mM=[0,5,10],
    library_strategy="RNA-Seq", library_selection="polyA", library_layout="PAIRED",
    instrument="Illumina NovaSeq (<TODO confirm>)", strandedness="unstranded",
    reference_genome="GRCh38", annotation="GENCODE v46",
    conditions="Mock/Cd x butyrate 0/5/10 mM, 3 replicates",
    provider="Novogene", vocab_notes="terms aligned to EBI ArrayExpress / MIABIS where possible",
  ),
  "DECIDE_A04_epithelial-coinfection": dict(
    title="Epithelial co-infection model (B. thetaiotaomicron + C. difficile), HT29",
    decide_project="A04",
    organism="Homo sapiens", ncbi_taxon="NCBI:txid9606",
    tissue="colon (epithelial cell line)", tissue_ontology="UBERON:0001155",
    model_system="HT29 colorectal adenocarcinoma epithelial cell line",
    cell_line="HT29", donor=None,
    pathogen="Clostridioides difficile", pathogen_taxon="NCBI:txid1496",
    commensal="Bacteroides thetaiotaomicron", commensal_taxon="NCBI:txid818",
    treatment=None, treatment_doses_mM=None,
    library_strategy="RNA-Seq", library_selection="polyA", library_layout="PAIRED",
    instrument="<TODO confirm>", strandedness="<TODO confirm>",
    reference_genome="GRCh38", annotation="GENCODE (<TODO version>)",
    conditions="Mock / Bt / Cd / Co (co-colonization)",
    provider="<TODO>", vocab_notes="terms aligned to EBI ArrayExpress / MIABIS where possible",
  ),
  "DECIDE_A03_MercedesGomez": dict(
    title="Molecular decision points of neonatal skin colonization and infection risk",
    decide_project="A03",
    organism="Mus musculus", ncbi_taxon="NCBI:txid10090",
    tissue="skin", tissue_ontology="UBERON:0002097",
    model_system="neonatal mouse (GF / SPF colonization)",
    cell_line=None, donor=None,
    pathogen=None, pathogen_taxon=None,
    treatment="bacterial/phage colonization", treatment_doses_mM=None,
    library_strategy="RNA-Seq", library_selection="polyA", library_layout="PAIRED",
    instrument="<TODO confirm>", strandedness="<TODO confirm>",
    reference_genome="GRCm39/mouse", annotation="GENCODE mouse (<TODO version>)",
    conditions="GF / SPF x bacteria / bacteria+phage",
    provider="<TODO>", vocab_notes="terms aligned to EBI ArrayExpress / MIABIS where possible",
  ),
}
cfg = STUDIES.get(NAME, dict(title=NAME, decide_project="<TODO>", organism="<TODO>"))

# ---------------- scan raw fastqs ----------------
patterns = ["**/*.fastq.gz","**/*.fq.gz","**/*.fastq","**/*.fq"]
files = []
for p in patterns:
    files += glob.glob(str(RAW / p), recursive=True)
files = sorted(set(files))

def sample_id(fn):
    b = os.path.basename(fn)
    # organoid: M0R1_R1.fastq.gz / Cd10R2_R2...  epithelial/other: <sample>_<...>_R1
    m = re.match(r'((?:M|Cd)\d+R\d+)', b)
    if m: return m.group(1)
    m = re.match(r'([A-Za-z0-9\-]+?)_(?:R?[12]|MKRN|S\d)', b)
    return m.group(1) if m else re.sub(r'(_R?[12])?\.(fastq|fq)\.gz$','',b)

def read_mate(fn):
    b=os.path.basename(fn)
    if re.search(r'_R?1[._]', b) or b.endswith('_1.fq.gz') or '_1.f' in b: return "R1"
    if re.search(r'_R?2[._]', b) or b.endswith('_2.fq.gz') or '_2.f' in b: return "R2"
    return "?"

from collections import defaultdict
samples = defaultdict(list)
for f in files:
    samples[sample_id(f)].append(f)

# ---------------- infer per-sample condition (organoid) ----------------
def cond_of(sid):
    m=re.match(r'(M|Cd)(\d+)R(\d+)', sid)
    if m:
        inf = "Mock" if m.group(1)=="M" else "Cd_infected"
        return dict(infection=inf, butyrate_mM=int(m.group(2)), replicate=int(m.group(3)))
    return dict(condition="<TODO>", replicate="<TODO>")

# ---------------- write YAML (hand-rolled, no yaml dep needed) ----------------
def y(v, ind=0):
    sp="  "*ind
    if isinstance(v, dict):
        out=[]
        for k,val in v.items():
            if isinstance(val,(dict,list)) and val:
                out.append(f"{sp}{k}:"); out.append(y(val, ind+1))
            else:
                out.append(f"{sp}{k}: {fmt(val)}")
        return "\n".join(out)
    if isinstance(v, list):
        return "\n".join(f"{sp}- {fmt(i)}" if not isinstance(i,dict) else f"{sp}-\n"+y(i,ind+1) for i in v)
    return f"{sp}{fmt(v)}"
def fmt(v):
    if v is None: return "null"
    if isinstance(v,bool): return str(v).lower()
    if isinstance(v,(int,float)): return str(v)
    s=str(v)
    return f'"{s}"' if (":" in s or s.strip()!=s or s=="") else s

header = {
  "fair_metadata_version": "0.1-draft",
  "generated": date.today().isoformat(),
  "note": "Draft FAIR metadata for review. <TODO> = confirm before publication.",
  "dataset": {
    "title": cfg.get("title"),
    "decide_project": cfg.get("decide_project"),
    "controlled_vocabulary": cfg.get("vocab_notes","EBI ArrayExpress / MIABIS aligned"),
    "biological": {k:cfg.get(k) for k in
        ["organism","ncbi_taxon","tissue","tissue_ontology","model_system","cell_line",
         "donor","pathogen","pathogen_taxon","commensal","commensal_taxon","treatment",
         "treatment_doses_mM","conditions"] if k in cfg},
    "technical": {k:cfg.get(k) for k in
        ["library_strategy","library_selection","library_layout","instrument","strandedness",
         "reference_genome","annotation","provider"] if k in cfg},
  }
}
sample_block=[]
for sid in sorted(samples):
    fq = sorted(samples[sid])
    rec = {"sample_id": sid}
    rec.update(cond_of(sid))
    rec["files"] = [{
        "path": os.path.relpath(f, ROOT),
        "read": read_mate(f),
        "metadata_ref": f"data/meta_data/{NAME}_raw_metadata.yaml#/samples/{sid}",
    } for f in fq]
    sample_block.append(rec)

out = y(header) + "\n\nsamples:\n" + y({r["sample_id"]:r for r in sample_block}, 1) + "\n"
outfile = META / f"{NAME}_raw_metadata.yaml"
outfile.write_text(out)
print(f"wrote {outfile}")
print(f"  samples: {len(samples)} | fastq files: {len(files)}")
if len(files)==0:
    print("  NOTE: no fastq found under data/raw_data (symlinks/tar only?) — file pointers will be empty.")
    print("        raw may be inside the tar; pointers still list expected sample IDs from names if present.")
