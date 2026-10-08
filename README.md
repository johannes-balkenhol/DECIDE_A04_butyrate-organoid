# DECIDE_A04 — C. difficile infection of colonic organoids under a butyrate dose-response

**Consortium:** DECIDE · Project Area A · sub-project A04
**Group:** Faber (collaborators: Sara, Alex, Franzi)
**Analyst:** Johannes Balkenhol
**Status:** active — downstream analysis + collaborator figures delivered; report drafted

---

## TL;DR — what this project is

Patient-derived human **colonic organoids** were infected with *Clostridioides difficile* (**Cd**) across a
**sodium-butyrate dose-response (0 / 5 / 10 mM)**, crossed with Mock. Design: Mock/Cd × butyrate 0/5/10 mM,
3 replicates = **18 libraries** (M0/M5/M10, Cd0/Cd5/Cd10 · R1–R3). Host mRNA mapped to GRCh38 (GENCODE v46);
differential expression with DESeq2/pydeseq2. Key readout: the **NF-κB / IL-8 (CXCL8) inflammatory programme**
and its coupling to an ER-chaperone arm. Main finding: butyrate acts as a **state-dependent brake** on IL-8 —
disproportionately effective during infection (dose-gated super-additive suppression at 10 mM), consistent
with an HDAC-inhibition mechanism on infection-opened NF-κB chromatin.

## Data lifecycle (FAIR)
`data/raw_data` (immutable FASTQ + Novogene delivery, kept on HPC) → `data/primary_data` (cleaned) →
`data/processed_data` (counts, DE) → `results/` (figures, tables). One-way, code-only.

## Structure
- `data/` — raw_data (FASTQ; large raw archived, on HPC), primary_data, processed_data, meta_data, refs, var
- `documents/` — SOPs, methods, DMP, QC report, reports, planning
- `scripts/` — pipeline steps (00_setup … 06_pathway) + `notebooks/` (numbered: QC, DE, dose-response,
  interaction, pathway, decision-points, IL-8 cross-model, collaborator figures)
- `analysis/` — per-step working outputs
- `results/` — figures, tables, DE_export_all, overview_for_sara

## Reproduce
```bash
conda env create -f environment.yml   # env: rnaseq (pydeseq2, scanpy, matplotlib)
# run scripts/notebooks/ in numerical order
```

## FAIR checklist
- [x] Canonical folder structure
- [x] README (this file) describing the project
- [x] DMP.md present at root
- [x] LICENSE at root  (open license — target: as permissive as consortium allows)
- [x] Raw-data metadata (`data/meta_data/DECIDE_A04_CdiffOrganoid_Butyrate_raw_metadata.yaml`), 18 samples / 82 FASTQ pointers
- [x] environment.yml for reproducibility
- [ ] Metadata harmonised to ArrayExpress + HubMap required fields  (in progress)
- [ ] Raw synced to cloud (currently HPC-only; primary_data can be synced for review)  <decision>
- [ ] Target repository accession (ArrayExpress)  <TODO>

## Data management
- **DMP:** see `DMP.md` (root).
- **Raw data location:** on HPC (`data/raw_data/`, incl. Novogene tar + extraction). Not synced to cloud by
  default (size). Primary/processed + code + results are synced for review.
- **Target publication repository:** ArrayExpress.
- **License:** open license (to confirm).
- **ELN / LabFolder:** <TODO add link>
- **Donor / ethics:** single donor <TODO donor ID / ethics ref>
