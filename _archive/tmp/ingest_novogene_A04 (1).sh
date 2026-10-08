#!/usr/bin/env bash
# =============================================================================
# ingest_novogene_A04.sh
# Unpack + verify + link a Novogene WinSCP delivery into data/raw/.
# Run from the project root:
#   cd ~/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate
#   bash ingest_novogene_A04.sh
# Safe to re-run. Does NOT delete the original tar.
# =============================================================================
set -euo pipefail

PROJ="$(pwd)"
RAW="${PROJ}/data/raw"
SAMPLES="${PROJ}/configs/samples.tsv"
mkdir -p "${RAW}" "${PROJ}/configs"

# ---- 0. guard: upload still running? ----------------------------------------
if ls "${PROJ}"/*.filepart >/dev/null 2>&1; then
  echo "Upload still in progress (.filepart present). Wait until it's gone, then re-run."
  exit 0
fi
TAR="$(ls "${PROJ}"/*.tar 2>/dev/null | head -1 || true)"
[[ -z "${TAR}" ]] && { echo "No .tar found in ${PROJ}. Nothing to ingest."; exit 1; }
echo ">> delivery tarball: ${TAR}"

# ---- 1. verify the tar (if a top-level checksum is provided) -----------------
# Novogene sometimes ships an MD5 for the tar itself; if MD5.txt lists the tar, check it.
if [[ -f "${PROJ}/MD5.txt" ]] && grep -q "$(basename "${TAR}")" "${PROJ}/MD5.txt"; then
  echo ">> verifying tarball against MD5.txt ..."
  ( cd "${PROJ}" && grep "$(basename "${TAR}")" MD5.txt | md5sum -c - ) \
    || { echo "TARBALL MD5 MISMATCH — re-upload, do not proceed."; exit 1; }
fi

# ---- 2. extract -------------------------------------------------------------
EXDIR="${PROJ}/_novogene_extracted"
mkdir -p "${EXDIR}"
echo ">> extracting (this can take a while for ~80 GB) ..."
tar -xf "${TAR}" -C "${EXDIR}"

# ---- 3. verify per-file checksums inside the delivery -----------------------
# Novogene puts an MD5.txt (often per-sample or one global) alongside the .fq.gz files.
echo ">> verifying per-file md5 (every MD5.txt found in the delivery) ..."
fail=0
while IFS= read -r md5file; do
  d="$(dirname "${md5file}")"
  ( cd "${d}" && md5sum -c "$(basename "${md5file}")" ) || fail=1
done < <(find "${EXDIR}" -name 'MD5.txt')
[[ "${fail}" -ne 0 ]] && { echo "ONE OR MORE MD5 CHECKS FAILED — investigate before analysis."; exit 1; }
echo "   all md5 checks passed."

# ---- 4. inventory FASTQ -----------------------------------------------------
mapfile -t R1 < <(find "${EXDIR}" -name '*_1.fq.gz' -o -name '*_R1*.f*q.gz' | sort)
echo ">> found ${#R1[@]} R1 files."
[[ ${#R1[@]} -eq 0 ]] && { echo "No FASTQ found under ${EXDIR}. Check the layout manually."; exit 1; }

# ---- 5. link into data/raw/ matching samples.tsv sample_ids -----------------
# Strategy: for each sample_id (M0R1..Cd10R3), find an R1 whose filename contains it.
if [[ ! -f "${SAMPLES}" ]]; then
  echo "WARN: ${SAMPLES} missing — dropping the populated samples.tsv there first is recommended."
  echo "      Listing FASTQ so you can map manually:"
  printf '   %s\n' "${R1[@]}"
  exit 0
fi

echo ">> linking FASTQ into data/raw/ by sample_id match ..."
unmatched=()
tail -n +2 "${SAMPLES}" | grep -v '^#' | cut -f1 | while read -r sid; do
  [[ -z "${sid}" ]] && continue
  r1="$(find "${EXDIR}" -name "*${sid}*_1.fq.gz" -o -name "*${sid}*_R1*.f*q.gz" 2>/dev/null | head -1)"
  r2="$(find "${EXDIR}" -name "*${sid}*_2.fq.gz" -o -name "*${sid}*_R2*.f*q.gz" 2>/dev/null | head -1)"
  if [[ -n "${r1}" && -n "${r2}" ]]; then
    ln -sf "${r1}" "${RAW}/${sid}_R1.fastq.gz"
    ln -sf "${r2}" "${RAW}/${sid}_R2.fastq.gz"
    echo "   linked ${sid}"
  else
    echo "   !! no match for ${sid}"
  fi
done

echo ""
echo ">> data/raw/ now contains:"
ls -l "${RAW}" | sed 's/^/   /'
echo ""
echo "If any sample shows '!! no match', the Novogene filenames don't contain the M/Cd IDs."
echo "In that case paste the output of:  find ${EXDIR} -name '*_1.fq.gz' | sort"
echo "and update the fastq_r1/fastq_r2 columns of samples.tsv to the real names."
