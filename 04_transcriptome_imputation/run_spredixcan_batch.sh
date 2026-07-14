#!/usr/bin/env bash
# 04_transcriptome_imputation/run_spredixcan_batch.sh
#
# Runs S-PrediXcan across all GTEx v8 MASHR tissue models for one or more
# recoded GWAS traits. Resumable: skips any (study, tissue) pair whose output
# file already exists, so an interrupted batch can simply be re-launched.
#
# Prerequisites:
#   - MetaXcan/software checked out, with the numpy2 compatibility patches
#     in patches/ applied to metax/metaxcan/AssociationCalculation.py
#   - GWAS tables recoded to varID via recode_gwas_to_varid.py
#   - GTEx v8 MASHR models downloaded (mashr_<tissue>.db + .txt.gz per tissue)
#
# Usage:
#   ./run_spredixcan_batch.sh <metaxcan_software_dir> <models_dir> <gwas_dir> <output_dir> <study1> [<study2> ...]

set -uo pipefail

SOFT="$1"; MODELS="$2"; GWASDIR="$3"; OUTDIR="$4"
shift 4
STUDIES=("$@")

mkdir -p "$OUTDIR"
cd "$SOFT"

n=0
total=$(( $(ls "$MODELS"/*.db | wc -l) * ${#STUDIES[@]} ))
for db in "$MODELS"/*.db; do
  tissue=$(basename "$db" .db | sed 's/^mashr_//')
  cov="$MODELS/mashr_${tissue}.txt.gz"
  for study in "${STUDIES[@]}"; do
    n=$((n+1))
    out="$OUTDIR/${study}__${tissue}.csv"
    if [[ -f "$out" ]]; then
      echo "[$n/$total] SKIP exists: $(basename "$out")"
      continue
    fi
    echo "[$n/$total] Running: study=$study tissue=$tissue"
    python3 SPrediXcan.py \
      --model_db_path "$db" \
      --model_db_snp_key varID \
      --covariance "$cov" \
      --gwas_file "$GWASDIR/${study}.recoded.table.gz" \
      --snp_column variant_id \
      --effect_allele_column effect_allele \
      --non_effect_allele_column non_effect_allele \
      --beta_column effect_size \
      --se_column standard_error \
      --pvalue_column pvalue \
      --keep_non_rsid \
      --output_file "$out" \
      --verbosity 6
  done
done
echo "BATCH_DONE"
