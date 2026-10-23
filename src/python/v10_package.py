#!/usr/bin/env python3
"""Package v10 primary audit outputs; no source mutation, no ranking by annotation."""
import os
import csv,json,hashlib,collections
from pathlib import Path
P=Path(os.environ.get("LOQUAT_PROJECT_ROOT",str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
O=P/"11_methods_audit_20261020"
V=O/"v10_universal_synteny"
def read(x):
 with open(x,newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
orig={x["apple_geneid"]:x for x in read(P/"07_loquat_projection/cross_species/all_strict_ortholog_scored.csv")}
rows=read(V/"s0.8_b0_external.csv")
pool=read(V/"s0.8_b0_pool.csv")
both=read(V/"s0.8_b0_both.csv")
assert len(pool)==190 and len(rows)==24 and len(both)==5,(len(pool),len(rows),len(both))
def is_true(x):return str(x).lower()=="true"
assert all(is_true(x["eligible_expression_direction"]) for x in pool)
assert all(is_true(x["external_any"]) for x in rows)
assert all(is_true(x["external_both"]) for x in both)
assert all(int(x["n_reference_apple_anchors"])>=8 and float(x["synteny_score_v10_min"])>=.8 for x in pool)
summary=json.loads((V/"summary.json").read_text())
annotation=[]
for r in rows:
 source=orig[r["apple_geneid"]]
 annotation.append({**r,
    "apple_gene_name_POST_SELECTION":source["apple_name"],
    "pear_description_POST_SELECTION":source["pear_description"],
    "functional_categories_POST_SELECTION":source["categories"],
    "apple_stage1_log2FC":source["apple_lfc_s1"],
    "apple_stage2_log2FC":source["apple_lfc_s2"],
    "original_legacy_score_NOT_USED":source["score"]})
annotation.sort(key=lambda r:(not is_true(r["external_both"]),-float(r["synteny_score_v10_min"]),r["apple_geneid"]))
with open(O/"V10_external24_annotated_AFTER_selection.csv","w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=annotation[0].keys());w.writeheader();w.writerows(annotation)
oldpool={r["apple_geneid"] for r in read(P/"09_tables/submission_v9/TableS1_candidate_pool38.csv")}
oldext={r["apple_geneid"] for r in read(P/"09_tables/submission_v9/TableS2_external_supported17.csv")}
newpool={r["apple_geneid"] for r in pool}
newext={r["apple_geneid"] for r in rows}
diff=[{"category":"pool_old_lost","gene":g} for g in sorted(oldpool-newpool)]+[{"category":"pool_new_added","gene":g} for g in sorted(newpool-oldpool)]+[{"category":"external_old_lost","gene":g} for g in sorted(oldext-newext)]+[{"category":"external_new_added","gene":g} for g in sorted(newext-oldext)]
with open(O/"V10_candidate_differences.csv","w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=["category","gene"]);w.writeheader();w.writerows(diff)
manifest=[
 P/"01_metadata/apple_main36_samples.tsv",
 P/"07_loquat_projection/cross_species/all_strict_ortholog_scored.csv",
 P/"09_tables/Table_candidate_blind_ranking_all20795.csv",
 P/"09_tables/submission_v9/TableS1_candidate_pool38.csv",
 P/"09_tables/submission_v9/TableS2_external_supported17.csv",
 P/"03_reference/apple_GDDH13/GCF_002114115.1_ASM211411v1_feature_table.txt.gz",
 P/"03_reference/pear_Cuiguan_GWHBAOS/GWHBAOS00000000.gff.gz",
 P/"03_reference/loquat_TBY/GWHETRZ00000000.1.gff.gz",
 O/"KEGG_frozen_TERM2GENE.csv",
 O/"KEGG_frozen_mapped_background.txt",
 O/"v10_universal_synteny.py",
 O/"KEGG_v10_threshold_sensitivity.R",
 O/"KEGG_v10_PKS_leaveout.R",
 O/"v10_package.py",
]
checks=[]
for f in manifest:
 # A processed-data release can omit original reference mirrors or local scripts.
 # Explicitly mark these as absent instead of aborting after computing v10 results.
 if not f.is_file():
  checks.append((str(f.relative_to(P)),"MISSING_IN_PARTIAL_PUBLIC_WORKSPACE",0))
  continue
 h=hashlib.sha256()
 with f.open("rb") as fp:
  while (b:=fp.read(1024*1024)):h.update(b)
 checks.append((str(f.relative_to(P)),h.hexdigest(),f.stat().st_size))
with open(O/"V10_SOURCE_SHA256.tsv","w") as f:
 f.write("source_path\tSHA256\tbytes\n")
 for row in checks:f.write("\t".join(map(str,row))+"\n")
fl=read(O/"KEGG_v10_flavonoid_summary.csv")
x={r["id"]:r for r in fl}
report={
 "pipeline":"v10 reproducible audit candidate universe",
 "primary_parameters":{"anchors_per_side":5,"target_window_strict_RBH_ranks":10,"collinearity_lis_orientation_free":True,
                       "combine_species":"minimum","minimum_available_apple_anchors":8,"primary_threshold":0.8,
                       "annotation_blind":True,"ranking_function_annotation_used":False},
 "counts":{"eligible":269,"pool":len(pool),"external_any":len(rows),"external_same_direction_both_tissues":len(both)},
 "flavonoid_fixed_151_pathway_FDR":{
  k:float(x[k]["FDR_BH_151"]) for k in ("s0.8_b0_pool","s0.8_b0_external","s0.8_b0_both")},
 "flavonoid_hits":{k:int(x[k]["hits"]) for k in ("s0.8_b0_pool","s0.8_b0_external","s0.8_b0_both")},
 "relative_to_old":{"pool_intersection":len(oldpool&newpool),"external_intersection":len(oldext&newext),
                    "old_external_lost":sorted(oldext-newext),"external_added":sorted(newext-oldext)},
 "limits":["Biological independence of extraction pools remains unclear","Revised synteny algorithm is a pre-specified exploratory audit, not a reproduction of missing legacy code",
           "KEGG membership is frozen from 2026-10-20 and differs from historical KEGG release",
           "Predicted RBH orthology and expression are not proof of causal dwarfing regulators"]}
leaveout=read(O/"KEGG_v10_PKS_leaveout_summary.csv")
report["flavonoid_PKS_leaveout_FDR"]={f"{r['group']}_{r['leaveout']}":float(r['FDR_BH_151']) for r in leaveout}
(O/"V10_RESULTS.json").write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps(report,indent=2,ensure_ascii=False))
