#!/usr/bin/env python3
"""Construct a *writable isolated* analysis workspace from unpacked Zenodo inputs.

No files in the dataset are overwritten by main DESeq2, v10 gene selection,
KEGG ORA, figure, or Supplementary scripts. A few original genome annotations
must still be separately obtained for full candidate score recomputation.
"""
from pathlib import Path
import argparse,json,shutil

def symlink(src, dest):
    src=Path(src).expanduser().resolve()
    if not src.exists():raise FileNotFoundError(str(src))
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists() or dest.is_symlink():raise FileExistsError(str(dest))
    dest.symlink_to(src,target_is_directory=src.is_dir())

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset",required=True,help="Unpacked Zenodo Loquat_Processed_Data_v1.0.0")
    ap.add_argument("--workspace",required=True,help="NEW EMPTY directory, separate from dataset")
    ap.add_argument("--apple-feature-table",help="Externally obtained GCF_002114115.1 feature_table.txt.gz")
    ap.add_argument("--pear-gff",help="Externally obtained GWHBAOS00000000.gff.gz")
    ap.add_argument("--loquat-gff",help="Externally obtained GWHETRZ00000000.1.gff.gz")
    ap.add_argument("--prefill-archived-results",action="store_true",
        help="Copy archived outputs to isolated workspace for figure-only regeneration; DO NOT treat as an independent DE rerun.")
    args=ap.parse_args()
    data=Path(args.dataset).expanduser().resolve()
    dest=Path(args.workspace).expanduser().resolve()
    if not (data/"DATA_DICTIONARY_HEADERS.tsv").is_file():
        raise ValueError("Dataset does not contain DATA_DICTIONARY_HEADERS.tsv")
    if dest==data or data in dest.parents or dest in data.parents:
        raise ValueError("Workspace must be outside of the source dataset")
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError("Refusing nonempty workspace; all outputs must be isolated")
    dest.mkdir(parents=True,exist_ok=True)
    mapping={
      "01_metadata":data/"01_sample_metadata",
      "04_quant":data/"02_salmon_quant",
      "12_transcript_to_gene_mapping":data/"12_transcript_to_gene_mapping",
      "07_KEGG_GSEA_full_gene_universe":data/"07_KEGG_GSEA_full_gene_universe",
      "07_loquat_projection/cross_species/all_strict_ortholog_scored.csv":data/"04_cross_species_mapping/all_strict_ortholog_scored.csv",
      "05_ortholog/rbh/apple_pear_loquat_RBH_triples.tsv":data/"04_cross_species_mapping/apple_centered_RBH_protein_triples.tsv",
      "09_tables/Table_candidate_blind_ranking_all20795.csv":data/"05_candidate_sets/Table_candidate_blind_ranking_all20795.csv",
      "09_tables/submission_v9/TableS1_candidate_pool38.csv":data/"11_historical_membership_for_audit_only/TableS1_candidate_pool38.csv",
      "09_tables/submission_v9/TableS2_external_supported17.csv":data/"11_historical_membership_for_audit_only/TableS2_external_supported17.csv",
      "12_submission_v10_assets/Figure3_Original_9Pathway_Source.csv":data/"10_figures/Figure3_Original_9Pathway_Source.csv",
      "12_submission_v10_assets/Figure4_24_heatmap_source.csv":data/"10_figures/Figure4_24_heatmap_source.csv",
      "12_submission_v10_assets/Table_S5_KEGG_fixed151_all_sets.csv":data/"10_figures/Table_S5_KEGG_fixed151_all_sets.csv",
      "12_submission_v10_assets/Table_S6_PKS_paralog_sensitivity.csv":data/"10_figures/Table_S6_PKS_paralog_sensitivity.csv",
      "12_submission_v10_assets/SuppFig1_V10_6_loci_anchors.json":data/"10_figures/SupplementaryFigureS1_6locus_anchors.json"}
    for name in ["KEGG_frozen_TERM2GENE.csv","KEGG_frozen_mapped_background.txt",
                 "KEGG_v10_synteny_full_FDR.csv","KEGG_v10_PKS_leaveout_ALL_paths.csv",
                 "KEGG_v10_PKS_leaveout_summary.csv","KEGG_v10_flavonoid_summary.csv","V10_RESULTS.json"]:
        mapping["11_methods_audit_20261020/"+name]=data/"06_KEGG_ORA_ortholog_restricted"/name
    ref={
      "03_reference/apple_GDDH13/GCF_002114115.1_ASM211411v1_feature_table.txt.gz":args.apple_feature_table,
      "03_reference/pear_Cuiguan_GWHBAOS/GWHBAOS00000000.gff.gz":args.pear_gff,
      "03_reference/loquat_TBY/GWHETRZ00000000.1.gff.gz":args.loquat_gff}
    for rel,p in ref.items():
        if p:mapping[rel]=Path(p)
    for rel,src in mapping.items():symlink(src,dest/rel)
    for rel in ["06_meta_analysis/apple_main36","06_meta_analysis/apple2026",
       "11_methods_audit_20261020/v10_universal_synteny",
       "12_submission_v10_assets","08_figures","03_reference"]:
        (dest/rel).mkdir(parents=True,exist_ok=True)
    if args.prefill_archived_results:
        for lab in ["apple_main36","apple2026"]:
            shutil.copytree(data/"03_differential_expression"/lab,
                            dest/"06_meta_analysis"/lab,dirs_exist_ok=True)
        for src in (data/"05_candidate_sets").iterdir():
            if src.is_file():
                shutil.copy2(src,dest/"11_methods_audit_20261020/v10_universal_synteny"/src.name)
    print(json.dumps({"workspace":str(dest),"prefilled_archived_results":args.prefill_archived_results,
     "linked_inputs":len(mapping),
     "raw_reference_annotations_complete":all(ref.values()),
     "safe_separate_output_directories":True},indent=2))
if __name__=="__main__":main()
