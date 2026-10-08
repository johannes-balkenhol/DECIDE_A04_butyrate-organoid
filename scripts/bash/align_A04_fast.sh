#!/usr/bin/env bash
# =============================================================================
# scripts/bash/align_A04_fast.sh
# 1) Rebuild data/raw/ correctly: per sample, MERGE all flowcell/lane FASTQ in
#    matched R1/R2 order (fixes the mismatched-mate bug from ingest's head -1).
# 2) Align all samples to GRCh38 in PARALLEL, tuned to ~100 cores.
#
# Lives in scripts/bash/ ; resolves the project root from its own location, so
# you can launch it from anywhere:
#   cd ~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate
#   nohup bash scripts/bash/align_A04_fast.sh > logs_align.out 2>&1 &
#   tail -f logs_align.out
#
# Tuning (defaults give 10 x 10 = 100 cores):
#   JOBS=10 THREADS=10  -> 100 cores   (default; 2 waves over 18 samples)
#   JOBS=12 THREADS=8   -> 96 cores
#   JOBS=7  THREADS=16  -> 112 cores   (fewer parallel, more per-sample)
# Note: STAR is CPU-only (no GPU); RAM ~30 GB per parallel job (10 -> ~300 GB,
# well within 2 TB). Past ~16 threads/sample STAR barely scales, so parallelism
# across samples is what fills the box.
# =============================================================================
set -euo pipefail

# ---- resolve project root from THIS script's location (scripts/bash/ -> ../..)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="${PROJ:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"

EXDIR="${PROJ}/_novogene_extracted"
RAW="${PROJ}/data/raw"
BAM="${PROJ}/data/processed/bam"
STAR_IDX="${PROJ}/data/refs/star_GRCh38"
SAMPLES="${PROJ}/configs/samples.tsv"
THREADS="${THREADS:-10}"      # threads per sample
JOBS="${JOBS:-10}"            # samples in parallel  (10 x 10 = 100 cores)

echo ">> project root: ${PROJ}"
[[ -f "${SAMPLES}" ]] || { echo "samples.tsv not found at ${SAMPLES} — wrong project root?"; exit 1; }
[[ -f "${STAR_IDX}/SAindex" ]] || { echo "STAR index missing at ${STAR_IDX} — run 'index' first."; exit 1; }
[[ -d "${EXDIR}" ]] || { echo "extracted data not found at ${EXDIR} — run ingest first."; exit 1; }
mkdir -p "${RAW}" "${BAM}"

log(){ echo -e "[$(date +%H:%M:%S)] $*"; }

# ---- 1. rebuild data/raw with correct per-sample lane merging ----------------
log "rebuilding data/raw/ with lane-aware merging ..."
find "${RAW}" -maxdepth 1 -type l -delete 2>/dev/null || true   # drop old (bad) symlinks; keep real merged files

tail -n +2 "${SAMPLES}" | grep -v '^#' | cut -f1 | while read -r sid; do
  [[ -z "${sid}" ]] && continue
  out1="${RAW}/${sid}_R1.fastq.gz"; out2="${RAW}/${sid}_R2.fastq.gz"
  if [[ -s "${out1}" && -s "${out2}" ]]; then echo "   keep ${sid} (already built)"; continue; fi

  sdir="$(find "${EXDIR}" -type d -name "${sid}" | head -1)"
  [[ -z "${sdir}" ]] && { echo "   !! no extracted dir for ${sid}"; continue; }

  mapfile -t r1 < <(ls "${sdir}"/*_1.fq.gz 2>/dev/null | sort)   # sort both the SAME way ->
  mapfile -t r2 < <(ls "${sdir}"/*_2.fq.gz 2>/dev/null | sort)   # R1/R2 lane order matches -> mates paired
  if [[ ${#r1[@]} -ne ${#r2[@]} || ${#r1[@]} -eq 0 ]]; then
    echo "   !! ${sid}: R1/R2 lane count mismatch (${#r1[@]} vs ${#r2[@]}) — check manually"; continue
  fi

  if [[ ${#r1[@]} -eq 1 ]]; then
    ln -sf "${r1[0]}" "${out1}"; ln -sf "${r2[0]}" "${out2}"
    echo "   ${sid}: 1 lane (symlink)"
  else
    echo "   ${sid}: ${#r1[@]} lanes -> concatenating (matched order)"
    cat "${r1[@]}" > "${out1}"
    cat "${r2[@]}" > "${out2}"
  fi
done

# ---- 2. parallel alignment (~100 cores) -------------------------------------
log "aligning: ${JOBS} samples in parallel x ${THREADS} threads = $((JOBS*THREADS)) cores ..."
export STAR_IDX BAM RAW THREADS
tail -n +2 "${SAMPLES}" | grep -v '^#' | cut -f1 | \
xargs -P "${JOBS}" -I{} bash -c '
  sid="$1"
  out="${BAM}/${sid}_"
  if [[ -s "${out}Aligned.sortedByCoord.out.bam" ]]; then echo "   skip ${sid} (done)"; exit 0; fi
  echo "   >> aligning ${sid}"
  STAR --runThreadN "${THREADS}" --genomeDir "${STAR_IDX}" \
    --readFilesIn "${RAW}/${sid}_R1.fastq.gz" "${RAW}/${sid}_R2.fastq.gz" \
    --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix "${out}" \
  && samtools index "${out}Aligned.sortedByCoord.out.bam" \
  && echo "   done ${sid}"
' _ {}

log "alignment finished. BAMs in ${BAM}/"
echo ">> next: bash run_pipeline_A04_hostonly.sh strand   # then quant"
