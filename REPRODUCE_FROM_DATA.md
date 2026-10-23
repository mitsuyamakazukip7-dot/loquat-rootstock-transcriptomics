# Reproduce from the companion Zenodo processed data

## 1. Create SAFE workspace (fresh directory outside the dataset)

python3 src/python/prepare_public_analysis_workspace.py --dataset /path/to/Loquat_Processed_Data_v1.0.0 --workspace /path/to/new_workspace --prefill-archived-results

This copies archived differential expression and candidate data into a separate writable output workspace for Figure2–5/S1 input compatibility. Pre-filling results is not a fresh DEG recomputation.

## 2. Independent quant.sf -> DEG recalculation

Use a different empty workspace without --prefill-archived-results, then:

LOQUAT_PROJECT_ROOT=/path/to/new_workspace Rscript src/R/apple_main36_deseq2.R

LOQUAT_PROJECT_ROOT=/path/to/new_workspace Rscript src/R/apple2026_deseq2.R

The main script uses a fixed full-gene KEGG table for the GSEA sensitivity rerun. The historical 37/27 manuscript panel originates from preserved historical data and the verified conditional-permutation workflow; it does not imply reselecting exactly 37 paths under arbitrary random seeds.

## 3. Recalculate local collinearity and candidate tiers

Obtain original GCF_002114115.1 NCBI feature_table.txt.gz, GWHBAOS00000000.gff.gz for pear, and GWHETRZ00000000.1.gff.gz for loquat (not mirrored in data archive).

Create another fresh workspace with the above command, including:
--apple-feature-table /path/to/apple_feature_table.txt.gz
--pear-gff /path/to/pear.gff.gz
--loquat-gff /path/to/loquat.gff.gz

Set LOQUAT_PROJECT_ROOT to this workspace, then execute:
python3 src/python/v10_universal_synteny.py
python3 src/python/build_supp_v10.py

Candidate selection: expression/direction eligible 269; strict local collinearity score >=0.80 AND 8 reference apple neighboring anchors minimum gives 190; 24 external-supported; five same-sign effects in both leaf contexts.

## 4. Pathway-permutation audit

cd verified_pathway_workflow && zsh RUN_FROZEN_AUDIT.sh

Expected: historical 37/27/0.395 conditional; frozen deterministic 38/28/0.395; all147 86/147/0.420.

## 5. Limits

Upstream reference acquisition/index-building, pear source-study full DEG table, DIAMOND and SRA Toolkit have additional requirements. A complete new-PC raw sequencing pipeline has NOT been tested. Historical scripts are supplied as research provenance, not a turnkey cross-platform installer.

Figure1 generated with the Python alternative is NOT the user's final BioRender conceptual Figure1. The figure license and final export remain pending. Data and code licenses remain author-approval pending.