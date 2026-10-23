#!/bin/zsh
set -euo pipefail
BASE="${0:A:h}"
export LOQUAT_PROJECT_ROOT="$BASE/input_snapshots"
export LOQUAT_AUDIT_DIR="$BASE"
cd "$BASE"
if [[ ! -s KEGG_MDM_FULL_TERM2GENE_20261022.csv ]]; then
 echo "ERROR: Frozen full KEGG mapping is absent. This package must run OFFLINE without a live KEGG query." >&2
 exit 3
fi
echo "Using project input root: $LOQUAT_PROJECT_ROOT"
echo "Using frozen KEGG: $BASE/KEGG_MDM_FULL_TERM2GENE_20261022.csv"
Rscript "$BASE/01_freeze_full_gsea_kegg.R" > "$BASE/recheck_step1.log" 2>&1
Rscript "$BASE/02_exact_permutation_full_frozen.R" > "$BASE/recheck_step2.log" 2>&1
Rscript "$BASE/03_joint_permutation_sensitivity.R" > "$BASE/recheck_step3.log" 2>&1
python3 "$BASE/VERIFY_RESULTS.py"
echo "ALL PASS. Optional seed sensitivity: Rscript 04_gsea_seed_sensitivity.R and 05_gsea_external_rng_seed_sensitivity.R"
