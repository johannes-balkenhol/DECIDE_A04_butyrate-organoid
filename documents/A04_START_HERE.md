# A04 — START HERE (single source of truth)

This chat is the one. Everything below reflects the real on-disk state on wbbi206 plus
the four prepared files. Ignore conflicting steps from any other chat.

## Locked facts (do not let these drift)
- **Host-only RNA-seq** (polyA). Align to **GRCh38 only**. *C. difficile* 630 (DSM27543)
  is methods text, NOT an alignment target. Do NOT download a Cd genome.
- **18 samples** = {Mock, Cd} x {0,5,10 mM butyrate} x 3 reps. Single donor (patient 35).
- **Model:** `~ infection_butyrate` group design. No donor term, no batch term.
- **Delivery:** Novogene tarball via WinSCP (not Nextcloud). Verify md5, extract, link.

## Current on-disk state (project root = ~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate)
- `X208SC26056228-Z01-F001.tar.filepart`  <- upload STILL RUNNING (wait for .filepart to vanish)
- `MD5.txt`, `checkSize.xlsx`, `QC Report/`  <- Novogene package (their own QC report)
- `tmp/organoid-cdiff-butyrate_sample data.xlsx`  <- the sample sheet
- `data/{raw, raw_data, processed, metadata, refs}`  <- note: raw_data is a stray dup, ignore it
- `scripts/{00_setup,01_qc,02_trim,03_align,04_quant,05_de,06_pathway}`  <- empty dirs from other chat; harmless
- empty: `notebooks/`, `results/`, `docs/`

## Place the 4 prepared files (from this chat)
| file | put it at |
|------|-----------|
| `samples.tsv` (18 rows, filled) | `configs/samples.tsv`  (run `mkdir -p configs` first) |
| `ingest_novogene_A04.sh` | project root |
| `run_pipeline_A04_hostonly.sh` | project root |
| `A04_downstream_DE.ipynb` | `notebooks/` |

(Optional: `A04_SOP_and_conventions_transfer.md` -> `docs/`.)

## Run order
**Now (data-independent, start immediately):**
```bash
cd ~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate
mkdir -p configs            # then drop samples.tsv into it
bash run_pipeline_A04_hostonly.sh index     # GRCh38+GENCODE download + STAR index (longest step)
```
**Read Novogene's QC while waiting** (saves running your own raw FastQC):
```bash
find "QC Report" -name '*_Report.html'      # open in browser; check Q30 (>85%), reads/sample, GC
head -5 MD5.txt ; ls -la *.tar* checkSize.xlsx
```
**When the upload finishes (.filepart gone):**
```bash
df -h /storage                              # ensure ~2x the tar size free (extraction doubles it)
bash ingest_novogene_A04.sh                 # md5 verify -> extract -> link FASTQ into data/raw/
bash run_pipeline_A04_hostonly.sh align     # STAR all 18 to GRCh38
bash run_pipeline_A04_hostonly.sh strand    # confirm strandedness from the data
bash run_pipeline_A04_hostonly.sh quant     # featureCounts -> results/tables/counts_raw.tsv
```
**Then:** open `notebooks/A04_downstream_DE.ipynb`, run top to bottom -> QC/PCA, DE per
contrast, volcanoes, A06 panel heatmap. Tables in `results/tables/DE_*.csv`.

## Skip these (other-chat steps that don't apply)
- `pull` step / Nextcloud download — data came via WinSCP.
- own raw FastQC first — Novogene's report covers it; only re-run if their summary looks off.
- `02_trim` — polyA usually needs none; trim only if QC shows adapter/quality problems.
- downloading the Cd genome — host-only.

## Loose ends to confirm (quick)
1. If `ingest` prints `!! no match` for samples, the Novogene filenames use library codes,
   not M0R1/Cd10R2. Paste `find _novogene_extracted -name '*_1.fq.gz' | sort` and map them.
2. `strand` step prints the strandedness call — set featureCounts `-s` to match.
3. After the transfer is done: tighten/remove `client_user` (it has group-write to the whole
   tree with a shared password).
4. Governance: this is patient-35 (pseudonymous) human data — keep the label, no re-identification.
