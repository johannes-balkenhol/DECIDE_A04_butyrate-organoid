#!/usr/bin/env bash
# =============================================================================
# verify_AR_absence.sh  — is AR truly unexpressed, or a reference/mapping artifact?
# Run from project root:  bash scripts/verify_AR_absence.sh
# Decides between: (1) genuine biology (AR off in colonic organoid)
#                  (2) AR missing/mis-annotated in THIS reference
#                  (3) reads present at locus but not counted (strand/multimap)
# =============================================================================
set -uo pipefail
PROJ="${HOME}/Projects_shared/DECIDE_A04_CdiffOrganoid_Butyrate"
GTF="${PROJ}/data/refs/human_gencode/gencode.v46.primary_assembly.annotation.gtf"
BAMDIR="${PROJ}/data/processed/bam"
AR_ID="ENSG00000169083"

echo "######## 1. Is AR ANNOTATED in the reference we built the index from? ########"
# If AR isn't in the GTF, STAR/featureCounts could never have counted it -> annotation artifact.
awk -F'\t' -v id="$AR_ID" '$3=="gene" && $9 ~ id {print "  GTF gene line:"; print "   "$1":"$4"-"$5" strand "$7; print "   "$9; found=1}
  END{if(!found) print "  !! AR gene_id "id" NOT in GTF -> annotation artifact, not biology"}' "$GTF"
echo "  exon count for AR in GTF:"
awk -F'\t' -v id="$AR_ID" '$3=="exon" && $9 ~ id {n++} END{print "   "n" exons annotated"}' "$GTF"

echo ""
echo "######## 2. Is the AR locus PHYSICALLY in the genome FASTA / index? ########"
FA="${PROJ}/data/refs/human_gencode/GRCh38.primary_assembly.genome.fa"
if [[ -f "${FA}.fai" || -f "${FA}" ]]; then
  [[ -f "${FA}.fai" ]] || samtools faidx "${FA}"
  echo "  chrX present in FASTA index?"; grep -P '^chrX\t' "${FA}.fai" | head -1 | awk '{print "   chrX length "$2}'
else
  echo "  (genome FASTA not found at ${FA} — skip; index was built though)"
fi

echo ""
echo "######## 3. Reads at the AR locus across ALL 18 BAMs (coverage + raw read count) ########"
# Pull AR coordinates from the GTF so there is no hard-coded-coords mismatch.
read CHR START END STRAND < <(awk -F'\t' -v id="$AR_ID" '$3=="gene" && $9 ~ id {print $1, $4, $5, $7; exit}' "$GTF")
echo "  AR locus from GTF: ${CHR}:${START}-${END} (${STRAND})"
# control: a gene we KNOW is on (TMPRSS2)
read C2 S2 E2 ST2 < <(awk -F'\t' '$3=="gene" && $9 ~ "gene_name \"TMPRSS2\"" {print $1,$4,$5,$7; exit}' "$GTF")
printf "  %-8s %12s %12s %12s\n" sample AR_reads AR_meandepth TMPRSS2_reads
for b in "${BAMDIR}"/*_Aligned.sortedByCoord.out.bam; do
  s=$(basename "$b" _Aligned.sortedByCoord.out.bam)
  arr=$(samtools view -c "$b" "${CHR}:${START}-${END}" 2>/dev/null)
  ard=$(samtools depth -r "${CHR}:${START}-${END}" "$b" 2>/dev/null | awk '{x+=$3;n++} END{printf "%.2f",(n?x/n:0)}')
  t2=$(samtools view -c "$b" "${C2}:${S2}-${E2}" 2>/dev/null)
  printf "  %-8s %12s %12s %12s\n" "$s" "$arr" "$ard" "$t2"
done

echo ""
echo "######## 4. Of any AR reads, are they being DISCARDED (multimap / wrong strand)? ########"
# If raw reads exist at the locus but counts=0, the cause is counting, not biology.
B1=$(ls "${BAMDIR}"/*_Aligned.sortedByCoord.out.bam | head -1)
echo "  on $(basename "$B1"):"
tot=$(samtools view -c "$B1" "${CHR}:${START}-${END}")
prim=$(samtools view -c -F 0x100 "$B1" "${CHR}:${START}-${END}")     # primary only
mapq=$(samtools view -c -q 255 "$B1" "${CHR}:${START}-${END}")       # STAR unique = MAPQ255
echo "   total reads at AR locus     : ${tot}"
echo "   primary alignments          : ${prim}"
echo "   STAR-unique (MAPQ 255)      : ${mapq}"
echo "  -> if total>0 but STAR-unique~0, AR reads are multimappers (counting issue)."
echo "  -> if total~0 everywhere, AR is genuinely not transcribed (biology). EXPECTED."

echo ""
echo "######## 5. Cross-check vs STAR's own per-gene counts (independent of featureCounts) ########"
RPG=$(ls "${BAMDIR}"/M0R1_ReadsPerGene.out.tab 2>/dev/null | head -1)
if [[ -n "${RPG:-}" ]]; then
  echo "  AR row in STAR ReadsPerGene (cols: unstr fwd rev) for M0R1:"
  grep "${AR_ID}" "$RPG" | head -1 | sed 's/^/   /' || echo "   AR not in STAR ReadsPerGene"
  echo "  (all three columns ~0 => not a featureCounts -s problem; AR simply has no reads)"
fi
echo ""
echo "DONE. Verdict logic:"
echo "  - AR in GTF + chrX in FASTA + locus reads~0 in all 18 + STAR counts~0  => GENUINE non-expression."
echo "  - AR not in GTF, or chrX absent                                        => reference artifact."
echo "  - locus reads>0 but unique~0 / STAR cols split                         => counting/strand artifact."
