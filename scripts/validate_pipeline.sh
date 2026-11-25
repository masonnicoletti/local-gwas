#!/usr/bin/env bash
# Compare pipeline outputs against the manual-run reference files.
set -u
P=results/pipeline
SUMMARY=$P/report/chr22_summary.tsv
CAUSAL=$P/report/chr22_causal.tsv
pass=0
fail=0

check_same() {
  if cmp -s "$1" "$2"; then
    echo "PASS  identical: $1"
    pass=$((pass + 1))
  else
    echo "FAIL  differs:   $1  vs  $2"
    fail=$((fail + 1))
  fi
}

check_value() {
  if [ "$2" = "$3" ]; then
    echo "PASS  $1 = $3"
    pass=$((pass + 1))
  else
    echo "FAIL  $1: got $2, expected $3"
    fail=$((fail + 1))
  fi
}

# value from the summary table: row label, column name
get() {
  awk -F'\t' -v r="$1" -v c="$2" \
    'NR==1 {for (i=1; i<=NF; i++) if ($i==c) col=i; next} $1==r {print $col}' "$SUMMARY"
}

for s in snps varqc sampleqc; do
  for ext in pgen pvar psam; do
    check_same "$P/plink/chr22_$s.$ext" "data/processed/chr22_$s.$ext"
  done
done

for k in in out; do
  check_same "$P/kinship/chr22_unrelated.king.cutoff.$k.id" "data/processed/chr22_unrelated.king.cutoff.$k.id"
  check_same "$P/pruning/chr22_pruned.prune.$k" "data/processed/chr22_pruned.prune.$k"
done

check_same "$P/pca/chr22_pca.eigenvec" results/pca/chr22_pca.eigenvec
check_same "$P/pca/chr22_pca.eigenval" results/pca/chr22_pca.eigenval

for m in nocov pcs; do
  check_same "$P/assoc/chr22_$m.trait.glm.linear" "results/assoc/chr22_$m.trait.glm.linear"
done

check_value "variants after conversion" "$(grep -vc '^#' $P/plink/chr22_snps.pvar)" 1055454
check_value "variants after variant QC" "$(grep -vc '^#' $P/plink/chr22_varqc.pvar)" 170949
check_value "samples after sample QC" "$(grep -vc '^#' $P/plink/chr22_sampleqc.psam)" 2422
check_value "variants kept by pruning" "$(wc -l < $P/pruning/chr22_pruned.prune.in | tr -d ' ')" 14744
check_value "lambda, no covariates" "$(get 'no covariates' lambda)" 2.041
check_value "lambda, with PCs" "$(get 'with PCs' lambda)" 1.076
check_value "significant, no covariates" "$(get 'no covariates' significant)" 423
check_value "significant, with PCs" "$(get 'with PCs' significant)" 81
check_value "causal detected, no covariates" "$(awk -F'\t' 'NR>1 && $1=="no covariates" && $NF=="True"' $CAUSAL | wc -l | tr -d ' ')" 6
check_value "causal detected, with PCs" "$(awk -F'\t' 'NR>1 && $1=="with PCs" && $NF=="True"' $CAUSAL | wc -l | tr -d ' ')" 6

echo
echo "passed: $pass   failed: $fail"
exit $((fail > 0))