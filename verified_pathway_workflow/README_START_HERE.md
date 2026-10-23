# LOQUAT KEGG GSEA permutation — self-contained frozen audit

**Status: verified local portable rerun on 2026-10-22.** No live KEGG requests are required.

## What is in this folder

- `input_snapshots/`: original biological-level VST expression matrix, primary apple DESeq2 stage results, two archived GSEA tables, the original 37-pathway selection table, and original permutation outputs. These are *inputs* from the frozen prior analysis.
- `KEGG_MDM_FULL_TERM2GENE_20261022.csv`: 21,003 genome-wide apple gene-pathway pairs across 163 pathways, retrieved 2026-10-22. Within the measured and size-filtered VST gene universe, 147 sets meet the 10–500 gene range.
- `01_freeze_full_gsea_kegg.R`: use archived KEGG mapping, fixed internal fgsea seed, recompute stage1/2 GSEA and compare with archived results. Mapping retrieval only happens if the supplied CSV is missing; `RUN_FROZEN_AUDIT.sh` explicitly refuses to run with it missing.
- `02_exact_permutation_full_frozen.R`: independently score 20 label assignments per stage, 400 combinations, using frozen 147 KEGG sets; compare historical37, new deterministic38, and all147 subsets.
- `03_joint_permutation_sensitivity.R`: sensitivity for a *hypothetical paired* stage design (20 coupled assignments) versus independent within-stage permutations (400).
- `04_gsea_seed_sensitivity.R` and `05_gsea_external_rng_seed_sensitivity.R`: stage1 borderline pathway (mdm00650) FDR sensitivity.
- `PERMUTATION_AUDIT_REPORT_JA.md`: interpretation, scope, caveats, and suggested manuscript wording.
- `VERIFY_RESULTS.py`: independent outcome assertions on the produced tables.

## One-command local replication

Requires R 4.5.1, clusterProfiler 4.18.4, limma 3.66.0, fgsea 1.36.2; Python3 for assertions, and no networking. Run in terminal after extracting:

```bash
zsh RUN_FROZEN_AUDIT.sh
```

A successful run prints `PASS portable frozen GSEA/permutation: 37/27/.395; 38/28/.395; 147/86/.420`.

## Main outputs

| Dataset | Reversal count | Tail fraction |
|---|---|---|
| Historical observed selected 37 | 27/37 | 158/400 = 0.395 |
| Deterministic full-GSEA selected 38 | 28/38 | 158/400 = 0.395 |
| All 147 fixed pathways | 86/147 | 168/400 = 0.420 |

The historical all147 tail was 170/400 = 0.425. The small change is attributable to pathway assignment/membership changes; 4 pathway sizes differ between archived GSEA and present full KEGG mapping.

**Statistical warning:** the selected 37 pathways are defined using the same original data. The 400-label comparison holds pathway identities fixed and does not reselect pathways in each permutation; therefore it is a *conditional descriptive sensitivity analysis*, **not** a selection-adjusted exact hypothesis test. The all147 analysis avoids observed-FDR pathway selection but still requires exchangeability of experimental units. Biological pooling and stage pairing in the original source study remain unresolved.

## Separate data/code archive

This folder reproduces **pathway-level GSEA and permutations from saved primary apple outputs**. The independent upstream 51-sample FASTQ→Salmon, DESeq2, RBH, candidate 269→190→24→5, and ORA reanalysis is documented separately in `../15_FULL_REPLICATION_20261021/` and `../14_DATA_CODE_RELEASE_AUDIT_20261021/`. The present audit package should not be misdescribed as a full 51-FASTQ archive.

No manuscript, original raw data, figure, or original analysis result was overwritten during the pathway audit.
