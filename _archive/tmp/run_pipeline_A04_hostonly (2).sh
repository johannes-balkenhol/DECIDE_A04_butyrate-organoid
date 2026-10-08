#!/usr/bin/env bash
# =============================================================================
# run_pipeline_A04_hostonly.sh
# DECIDE A04 — host-only RNA-seq (human colon organoid response to C. difficile,
# butyrate dose-response). polyA mRNA => align to GRCh38 ONLY. No bacterial lane.
#
# Run on wbbi206 inside the rnaseq env, from the project root:
#   cd ~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate
#   bash run_pipeline_A04_hostonly.sh <step>
# where <step> is one of: pull | qc | index | align | quant | strand | all
# Steps are safe to re-run; each skips work already done where practical.
# =============================================================================
set -euo pipefail

# ---- config -----------------------------------------------------------------
PROJ="${HOME}/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate"
THREADS="${THREADS:-16}"
GENCODE_RELEASE="${GENCODE_RELEASE:-46}"          # check gencodegenes.org for current
READ_LEN="${READ_LEN:-150}"                       # set from FASTQ; sjdbOverhang = READ_LEN-1
STRAND="${STRAND:-2}"                             # featureCounts: 0 unstr,1 fwd,2 rev(polyA TruSeq). VERIFY via 'strand' step.
NC_RAW="nextcloud:/DECIDE/Project Area A/A04/C-difficile-infection-of-organoid-model with-butyrate"

REFDIR="${PROJ}/data/refs/human_gencode"
GENOME_FA="${REFDIR}/GRCh38.primary_assembly.genome.fa"
GTF="${REFDIR}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf"
STAR_IDX="${PROJ}/data/refs/star_GRCh38"
RAW="${PROJ}/data/raw"
BAM="${PROJ}/data/processed/bam"
QC="${PROJ}/results/qc"
TAB="${PROJ}/results/tables"
SAMPLES="${PROJ}/configs/samples.tsv"

mkdir -p "${REFDIR}" "${STAR_IDX}" "${RAW}" "${BAM}" "${QC}" "${TAB}"

log(){ echo -e "\n[$(date +%H:%M:%S)] $*"; }

# ---- step: pull raw from consortium + verify --------------------------------
step_pull(){
  log "pulling raw from Nextcloud (host-only, ~86 GB) ..."
  rclone copy "${NC_RAW}/data/" "${RAW}/_incoming/" -P --transfers 4
  log "raw landed under ${RAW}/_incoming/ . Listing:"
  find "${RAW}/_incoming" -name '*.f*q.gz' | sort
  if [[ -f "${RAW}/_incoming/md5sums.txt" ]]; then
    log "verifying md5 ..."
    ( cd "${RAW}/_incoming" && md5sum -c md5sums.txt ) || { echo "MD5 MISMATCH — stop."; exit 1; }
  else
    echo "WARN: no md5sums.txt found — ask Sara for checksums before trusting the upload."
  fi
  echo ">> When verified, move/symlink the FASTQ into ${RAW}/ with names matching samples.tsv."
  echo "   (Novogene often names them <sample>_1.fq.gz / _2.fq.gz — adjust samples.tsv fastq cols.)"
}

# ---- step: raw QC -----------------------------------------------------------
step_qc(){
  log "FastQC on raw reads ..."
  mkdir -p "${QC}/fastqc_raw"
  fastqc -t "${THREADS}" -o "${QC}/fastqc_raw" "${RAW}"/*.f*q.gz
  multiqc -f -o "${QC}" "${QC}/fastqc_raw"
  log "MultiQC report: ${QC}/multiqc_report.html"
}

# ---- step: references + STAR index -----------------------------------------
step_index(){
  if [[ ! -s "${GENOME_FA}" ]]; then
    log "downloading GRCh38 + GENCODE v${GENCODE_RELEASE} ..."
    ( cd "${REFDIR}"
      wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/GRCh38.primary_assembly.genome.fa.gz"
      wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
      gunzip -kf GRCh38.primary_assembly.genome.fa.gz gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz )
  fi
  if [[ -f "${STAR_IDX}/SAindex" ]]; then
    log "STAR index already present — skipping build."
  else
    log "building STAR index (GRCh38, host-only) ..."
    STAR --runMode genomeGenerate --runThreadN "${THREADS}" \
      --genomeDir "${STAR_IDX}" \
      --genomeFastaFiles "${GENOME_FA}" \
      --sjdbGTFfile "${GTF}" \
      --sjdbOverhang $((READ_LEN-1))
  fi
}

# ---- step: align all samples ------------------------------------------------
step_align(){
  log "STAR alignment of 18 samples to GRCh38 ..."
  tail -n +2 "${SAMPLES}" | grep -v '^#' | while IFS=$'\t' read -r sid _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ r1 r2; do
    [[ -z "${sid:-}" ]] && continue
    local out="${BAM}/${sid}_"
    if [[ -s "${out}Aligned.sortedByCoord.out.bam" ]]; then
      echo "   skip (done): ${sid}"; continue
    fi
    echo "   aligning ${sid}"
    STAR --runThreadN "${THREADS}" --genomeDir "${STAR_IDX}" \
      --readFilesIn "${RAW}/${r1}" "${RAW}/${r2}" --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix "${out}"
    samtools index "${out}Aligned.sortedByCoord.out.bam"
  done
  multiqc -f -o "${QC}" "${BAM}"
}

# ---- step: verify strandedness (run after at least one BAM exists) ----------
step_strand(){
  log "checking strandedness from STAR ReadsPerGene (col4=fwd, col3=rev) ..."
  local f; f=$(ls "${BAM}"/*ReadsPerGene.out.tab | head -1)
  [[ -z "${f}" ]] && { echo "no ReadsPerGene yet — run align first."; exit 1; }
  awk 'NR>4{u+=$2; f+=$3; r+=$4} END{
    printf "  unstranded=%d  forward(col3)=%d  reverse(col4)=%d\n", u, f, r;
    if (r>2*f) print "  => REVERSE stranded: set STRAND=2 (featureCounts -s 2)";
    else if (f>2*r) print "  => FORWARD stranded: set STRAND=1";
    else print "  => UNSTRANDED: set STRAND=0" }' "${f}"
  echo "  (polyA TruSeq is usually reverse / -s 2; confirm against the printout above.)"
}

# ---- step: featureCounts ----------------------------------------------------
step_quant(){
  log "featureCounts (paired, -s ${STRAND}) ..."
  featureCounts -T "${THREADS}" -p --countReadPairs -s "${STRAND}" \
    -a "${GTF}" -o "${TAB}/counts_raw.tsv" \
    "${BAM}"/*_Aligned.sortedByCoord.out.bam
  log "counts: ${TAB}/counts_raw.tsv  (summary: ${TAB}/counts_raw.tsv.summary)"
  echo ">> sanity: assignment rate should be high; if very low, strandedness (-s) is wrong."
}

case "${1:-all}" in
  pull)   step_pull ;;
  qc)     step_qc ;;
  index)  step_index ;;
  align)  step_align ;;
  strand) step_strand ;;
  quant)  step_quant ;;
  all)    step_index; step_qc; step_align; step_strand; step_quant ;;
  *) echo "usage: $0 {pull|qc|index|align|strand|quant|all}"; exit 1 ;;
esac
log "done: ${1:-all}"
