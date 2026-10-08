#!/usr/bin/env bash
# ============================================================================
# DECIDE project harmonization — FAIR structure (bash / terminal)
# Conforms projects to the canonical A03 layout. DRY-RUN by default.
#
#   ./harmonize_decide.sh              # plan only, moves nothing (safe)
#   APPLY=1 ./harmonize_decide.sh      # actually move/mkdir
#   PROJECT=a04 APPLY=1 ./harmonize_decide.sh   # one project only
#
# Canonical layout (no src/, no code/, no documents/ — docs/ only;
# notebooks ONLY under scripts/notebooks/; dupes+unused -> _archive/):
#   data/{raw,processed,refs,metadata,var}
#   scripts/{notebooks,bash,00_setup,01_qc,02_trim,03_align,04_quant,05_de,06_pathway}
#   docs/  results/{figures,tables}  analysis/  figures/  archive/  _archive/
#   README.md environment.yml config.yaml CHANGELOG.md LICENSE
# ============================================================================
set -uo pipefail

# ---------------- CONFIG (edit these) ----------------
BASE="${PROJECTS_ROOT:-$HOME/Projects_shared}"
APPLY="${APPLY:-0}"                 # 0 = dry run, 1 = execute
ONLY="${PROJECT:-all}"              # all | a03 | a04 | epi

# project folder names on the HPC (source) and their harmonized names (target)
A03_SRC="DECIDE_A03_MercedesGomez";                    A03_DST="DECIDE_A03_neonatal-skin-colonization"
A04_SRC="DECIDE_A04_CdiffOrganoid_Butyrate";           A04_DST="DECIDE_A04_organoid-cdiff-butyrate"
EPI_SRC="epithelial_coinfection-model_WORKING";        EPI_DST="DECIDE_A04_epithelial-coinfection"
# NOTE: set A03_DST=$A03_SRC if you do NOT want to rename the reference folder.

CANON_DIRS=(
  data data/raw data/processed data/refs data/metadata data/var
  scripts scripts/notebooks scripts/bash
  scripts/00_setup scripts/01_qc scripts/02_trim scripts/03_align scripts/04_quant scripts/05_de scripts/06_pathway
  docs results results/figures results/tables analysis figures archive _archive
)

# ---------------- helpers ----------------
say(){ printf '%s\n' "$*"; }
run(){ # echo (dry) or execute
  if [[ "$APPLY" == "1" ]]; then eval "$@"; else say "  [dry] $*"; fi
}
ensure_dirs(){ local root="$1"; for d in "${CANON_DIRS[@]}"; do run "mkdir -p '$root/$d'"; done; }

# move a path into target slot, but if it already exists at destination, send to _archive
mv_safe(){ # mv_safe <src> <dst_dir> <root>
  local src="$1" dstdir="$2" root="$3"
  [[ -e "$src" ]] || return 0
  local name; name="$(basename "$src")"
  local dest="$root/$dstdir/$name"
  if [[ -e "$dest" ]]; then
    run "mkdir -p '$root/_archive/$dstdir'"
    run "mv -n '$src' '$root/_archive/$dstdir/${name}.dup'"
    say "  ~ dup -> _archive/$dstdir/${name}.dup"
  else
    run "mkdir -p '$root/$dstdir'"
    run "mv -n '$src' '$dest'"
  fi
}
archive(){ # archive <path> <root> [reason]
  local p="$1" root="$2" why="${3:-unused}"
  [[ -e "$p" ]] || return 0
  run "mkdir -p '$root/_archive'"
  run "mv -n '$p' '$root/_archive/'"
  say "  \u2192 _archive ($why): $(basename "$p")"
}

hr(){ say "------------------------------------------------------------"; }
header(){ hr; say "PROJECT: $1"; say "  src: $2"; say "  dst: $3   (APPLY=$APPLY)"; hr; }

# =====================================================================
# A03 — reference layout: already canonical. Only rename + minor tidy.
# =====================================================================
harmonize_a03(){
  local src="$BASE/$A03_SRC" dst="$BASE/$A03_DST"
  [[ -d "$src" ]] || { say "SKIP A03: $src not found"; return; }
  header "A03 (reference)" "$src" "$dst"
  if [[ "$A03_DST" != "$A03_SRC" ]]; then run "mv -n '$src' '$dst'"; else dst="$src"; fi
  ensure_dirs "$dst"
  # A03 is already clean; nothing else to move. (notebooks/ already canonical under top-level;
  # if you want them under scripts/notebooks, uncomment:)
  # [[ -d "$dst/notebooks" ]] && mv_safe "$dst/notebooks" "scripts" "$dst"
  say "A03 done (structure ensured)."
}

# =====================================================================
# A04 organoid — tidy: dup notebooks in tmp/, scripts/notebooks canonical,
# QC Report (space), data/raw_data dup, big raw handled by sync (not here).
# =====================================================================
harmonize_a04(){
  local src="$BASE/$A04_SRC" dst="$BASE/$A04_DST"
  [[ -d "$src" ]] || { say "SKIP A04: $src not found"; return; }
  header "A04 organoid" "$src" "$dst"
  if [[ "$A04_DST" != "$A04_SRC" ]]; then run "mv -n '$src' '$dst'"; else dst="$src"; fi
  ensure_dirs "$dst"

  # docs: A04 already has docs/ -> keep. configs/ -> keep at root (canonical has config.yaml; keep configs/ too)
  # tmp/ = scratch & duplicates -> _archive
  archive "$dst/tmp" "$dst" "scratch/duplicates"
  # data/raw_data empty duplicate of data/raw -> _archive
  [[ -d "$dst/data/raw_data" ]] && archive "$dst/data/raw_data" "$dst" "empty dup of data/raw"
  # 'QC Report' (space in name) -> docs/qc_report
  if [[ -d "$dst/QC Report" ]]; then run "mkdir -p '$dst/docs'"; run "mv -n '$dst/QC Report' '$dst/docs/qc_report'"; say "  QC Report -> docs/qc_report"; fi
  # loose top-level scripts -> scripts/ ; pipeline runner stays discoverable
  for f in ingest_novogene_A04.sh run_pipeline_A04_hostonly.sh align_A04_fast.sh; do
    [[ -e "$dst/$f" ]] && mv_safe "$dst/$f" "scripts/00_setup" "$dst"
  done
  # logs -> data/var
  [[ -e "$dst/logs_align.out" ]] && mv_safe "$dst/logs_align.out" "data/var" "$dst"
  # BIG RAW: leave the .tar and _novogene_extracted in place; sync excludes them (already in cloud).
  say "  (87GB tar + _novogene_extracted left in place; sync will EXCLUDE them — already in cloud)"
  say "A04 done."
}

# =====================================================================
# epithelial (ex-A06) — heavy dedup: ~40 causal_upstream_v* copies,
# case-colliding Networks/networks, parallel scripts/{exploratory,notebooks,core_pipeline}.
# Keep ONE canonical notebook set under scripts/notebooks/, rest -> _archive.
# =====================================================================
harmonize_epi(){
  local src="$BASE/$EPI_SRC" dst="$BASE/$EPI_DST"
  [[ -d "$src" ]] || { say "SKIP epi: $src not found"; return; }
  header "epithelial co-infection (was A06 -> A04)" "$src" "$dst"
  if [[ "$EPI_DST" != "$EPI_SRC" ]]; then run "mv -n '$src' '$dst'"; else dst="$src"; fi
  ensure_dirs "$dst"

  # documents/ (non-canonical) -> docs/
  [[ -d "$dst/documents" ]] && { for f in "$dst"/documents/*; do mv_safe "$f" "docs" "$dst"; done; run "rmdir --ignore-fail-on-non-empty '$dst/documents'"; }

  # scripts: canonical set = scripts/core_pipeline (the curated one). Move it to scripts/notebooks.
  if [[ -d "$dst/scripts/core_pipeline" ]]; then
    for f in "$dst"/scripts/core_pipeline/*.ipynb "$dst"/scripts/core_pipeline/*.py; do
      [[ -e "$f" ]] && mv_safe "$f" "scripts/notebooks" "$dst"
    done
  fi
  # the parallel/duplicate notebook trees -> _archive (exploratory, notebooks, core_notebooks_*)
  for d in scripts/exploratory scripts/notebooks_old scripts/core_notebooks_20260608 scripts/notebooks/old_07012026; do
    [[ -e "$dst/$d" ]] && archive "$dst/$d" "$dst" "duplicate notebook tree"
  done
  # if scripts/notebooks already existed as a second parallel tree with v* dupes, archive the v* copies
  if [[ -d "$dst/scripts/notebooks" ]]; then
    for f in "$dst"/scripts/notebooks/10_causal_upstream_v*.ipynb "$dst"/scripts/notebooks/*"(1)"* "$dst"/scripts/notebooks/*"(2)"*; do
      [[ -e "$f" ]] && archive "$f" "$dst" "old version/duplicate"
    done
  fi
  # case-colliding Networks/ vs networks/ (dangerous on case-insensitive cloud): keep networks/, archive Networks/
  if [[ -d "$dst/analysis/Networks" && -d "$dst/analysis/networks" ]]; then
    archive "$dst/analysis/Networks" "$dst" "case collision with networks/"
  fi
  # same for GSEA/ vs gsea/
  if [[ -d "$dst/analysis/GSEA" && -d "$dst/analysis/gsea" ]]; then
    archive "$dst/analysis/GSEA" "$dst" "case collision with gsea/"
  fi
  # big loose zip at root -> _archive
  [[ -e "$dst/epithel_coinfection_analysis_20260608.zip" ]] && archive "$dst/epithel_coinfection_analysis_20260608.zip" "$dst" "large snapshot zip"
  [[ -e "$dst/untitled.txt" ]] && archive "$dst/untitled.txt" "$dst" "junk"
  # var/ (scratch) already exists -> keep as data/var
  [[ -d "$dst/var" ]] && { for f in "$dst"/var/*; do mv_safe "$f" "data/var" "$dst"; done; run "rmdir --ignore-fail-on-non-empty '$dst/var'"; }
  say "epithelial done."
}

# ---------------- run ----------------
say "DECIDE harmonization  |  BASE=$BASE  |  APPLY=$APPLY  |  ONLY=$ONLY"
say "(dry-run prints planned actions; set APPLY=1 to execute)"; hr
case "$ONLY" in
  a03) harmonize_a03 ;;
  a04) harmonize_a04 ;;
  epi) harmonize_epi ;;
  all) harmonize_a03; harmonize_a04; harmonize_epi ;;
  *) say "unknown PROJECT=$ONLY (use a03|a04|epi|all)";;
esac
hr; say "DONE. Review actions above."
[[ "$APPLY" != "1" ]] && say "This was a DRY RUN. Re-run with:  APPLY=1 ./harmonize_decide.sh"
