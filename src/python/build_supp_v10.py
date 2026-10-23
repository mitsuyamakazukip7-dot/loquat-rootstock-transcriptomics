#!/usr/bin/env python3
"""Submission-facing Supplementary CSVs from frozen v10 audit (no stale gene sets).
Build one source of truth for paper figures; all output in 12_submission_v10_assets.
"""
import csv, math, json, re
from pathlib import Path
import os
from collections import Counter
R=Path(os.environ.get("LOQUAT_PROJECT_ROOT",str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
A=R/"11_methods_audit_20261020"
V=A/"v10_universal_synteny"
O=R/"12_submission_v10_assets"
O.mkdir(parents=True,exist_ok=True)
def read(p):
 with open(p,newline="",encoding="utf-8-sig") as f: return list(csv.DictReader(f))
def save(name, data, keys=None):
 if not data:raise ValueError("Cannot save empty "+name)
 keys=keys or list(data[0])
 with open(O/name,"w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore");w.writeheader();w.writerows(data)
def number(v):
 try:
  x=float(v)
  return x if math.isfinite(x) else None
 except:return None
def signif(v):return (number(v) or 99)<.05
orig={str(r["apple_geneid"]):r for r in read(R/"07_loquat_projection/cross_species/all_strict_ortholog_scored.csv")}
new_de={r["GENEID"]:r for r in read(R/"06_meta_analysis/apple2026/DE_new_shoot_PYGD076_vs_M9T337.csv")}
spur_de={r["GENEID"]:r for r in read(R/"06_meta_analysis/apple2026/DE_fruiting_spur_PYGD076_vs_M9T337.csv")}
e=read(V/"eligible_269.csv")
keep={r["apple_geneid"] for r in read(V/"s0.8_b0_pool.csv")}
ext={r["apple_geneid"] for r in read(V/"s0.8_b0_external.csv")}
both={r["apple_geneid"] for r in read(V/"s0.8_b0_both.csv")}
legacy38={r["apple_geneid"] for r in read(R/"09_tables/submission_v9/TableS1_candidate_pool38.csv")}
legacy17={r["apple_geneid"] for r in read(R/"09_tables/submission_v9/TableS2_external_supported17.csv")}
assert len(e)==269 and len(keep)==190 and len(ext)==24 and len(both)==5
assert both<=ext<=keep
cols=[
"Apple_GeneID","Apple_Annotation_postselection","Pear_GeneID","Loquat_GeneID",
"Stage_Response_Class","Apple_Stage1_log2FC","Apple_Stage1_FDR",
"Apple_Stage2_log2FC","Apple_Stage2_FDR","Apple_Stage_x_Rootstock_FDR",
"Pear_log2FC","Pear_published_qvalue","Loquat_Stem_median_TPM",
"RBH_reference_anchors","Pear_collinear_anchors","Loquat_collinear_anchors",
"Pear_collinearity_fraction","Loquat_collinearity_fraction",
"V10_min_collinearity_score","Pass_V10_synteny_0p80",
"External_shoot_log2FC","External_shoot_FDR","External_spur_log2FC","External_spur_FDR",
"External_any_FDR_lt_0p05","External_both_same_direction_FDR_lt_0p05",
"V10_external24_member","V10_both5_member",
"Legacy38_member","Legacy17_member"]
out=[]
for q in e:
 r=orig[q["apple_geneid"]];g=q["apple_geneid"]
 d=dict(zip(cols,[
 g,r["apple_name"],q["pear_evm"],q["loquat_eja"],q["class"],
 number(r["apple_lfc_s1"]),number(r["apple_padj_s1"]),
 number(r["apple_lfc_s2"]),number(r["apple_padj_s2"]),number(r["apple_interaction_padj"]),
 number(r["pear_lfc"]),number(r["pear_qval"]),number(r["loq_tpm_median"]),
 int(q["n_reference_apple_anchors"]),int(q["pear_collinear_n"]),int(q["loquat_collinear_n"]),
 number(q["pear_collinear_fraction"]),number(q["loquat_collinear_fraction"]),
 number(q["synteny_score_v10_min"]),g in keep,
 number(new_de.get(g,{}).get("log2FoldChange")),number(r["validation_new_padj"]),
 number(spur_de.get(g,{}).get("log2FoldChange")),number(r["validation_fruit_padj"]),
 (signif(r["validation_new_padj"]) or signif(r["validation_fruit_padj"])),
 str(r["validation_both_same_dir"]).lower()=="true",
 g in ext,g in both,g in legacy38,g in legacy17]))
 out.append(d)
order={"stage_consistent":0,"interaction_associated":1}
out.sort(key=lambda x:(order[x["Stage_Response_Class"]],-x["V10_min_collinearity_score"],x["Apple_GeneID"]))
save("Table_S1_All_269_eligible.csv",out,cols)
save("Table_S2_Synteny_190.csv",[r for r in out if r["Pass_V10_synteny_0p80"]],cols)
save("Table_S3_External_24.csv",[r for r in out if r["V10_external24_member"]],cols)
save("Table_S4_BothTissues_5.csv",[r for r in out if r["V10_both5_member"]],cols)
# S5 uses ALL 151 pathways x 3 groups; no cherry-picking or zero-hit exclusion.
k=read(A/"KEGG_v10_synteny_full_FDR.csv")
tags={"s0.8_b0_pool":"S2_190","s0.8_b0_external":"S3_24","s0.8_b0_both":"S4_5"}
k2=[]
for t in k:
 if t["id"] not in tags:continue
 k2.append({"Candidate_Set":tags[t["id"]],"KEGG_ID":t["pathway"],
"Candidates_Total":int(float(t["original_n"])),"Candidates_KEGG_mapped":int(float(t["mapped_n"])),
"Genes_in_Pathway":int(float(t["hits"])),"Fold_enrichment":number(t["fold"]),
"Hypergeometric_p":number(t["raw_p"]),"BH_FDR_151_tested":number(t["FDR_BH_151"])})
k2.sort(key=lambda x:(x["Candidate_Set"],x["BH_FDR_151_tested"],x["KEGG_ID"]))
assert len(k2)==453
save("Table_S5_KEGG_fixed151_all_sets.csv",k2)
# PKS leave-out, full controls
p=read(A/"KEGG_v10_PKS_leaveout_summary.csv")
p2=[{"Candidate_Set":{"pool":"S2_190","external":"S3_24","both":"S4_5"}[r["group"]],
"PKS5_Exclusion":r["leaveout"],"Genes_in_Candidate_Set":int(r["raw_set_size"]),
"KEGG_mapped_genes":int(r["n_mapped"]),"Flavonoid_Gene_Hits":int(r["nHit"]),
"Raw_hypergeometric_p":number(r["p"]),"BH_FDR_151_tested":number(r["FDR_BH_151"])} for r in p]
save("Table_S6_PKS_paralog_sensitivity.csv",p2)
# Legacy membership audit: all genes from union old/new
changed=[r for r in out if r["Apple_GeneID"] in legacy38 or r["Pass_V10_synteny_0p80"]]
s7=[{"Apple_GeneID":r["Apple_GeneID"],"Apple_Annotation_postselection":r["Apple_Annotation_postselection"],
"Loquat_GeneID":r["Loquat_GeneID"],"V10_score":r["V10_min_collinearity_score"],
"Legacy38":r["Legacy38_member"],"V10_190":r["Pass_V10_synteny_0p80"],
"Legacy17":r["Legacy17_member"],"V10_24":r["V10_external24_member"],
"V10_both5":r["V10_both5_member"],
"Change":("retained" if r["Legacy38_member"] and r["Pass_V10_synteny_0p80"] else "new_in_v10" if r["Pass_V10_synteny_0p80"] else "excluded_by_v10")} for r in changed]
save("Table_S7_Legacy_vs_V10.csv",s7)
# Figure data from source of truth
fig=[r for r in out if r["V10_external24_member"]]
fig.sort(key=lambda x:(not x["V10_both5_member"],x["Stage_Response_Class"]!="stage_consistent",-x["V10_min_collinearity_score"],x["Apple_GeneID"]))
save("Figure4_24_heatmap_source.csv",fig)
fig5=[r for r in k2 if r["KEGG_ID"]=="mdm00941"]
save("Figure5_flavonoid_source.csv",fig5)
summary={
"strict_RBH_universe":20795,"eligible":len(out),"stage_consistent":sum(r["Stage_Response_Class"]=="stage_consistent" for r in out),"interaction_associated":sum(r["Stage_Response_Class"]=="interaction_associated" for r in out),
"collinearity_0p8_with_min8_anchors":len(keep),"external_any":len(ext),"external_two_tissues_same_direction":len(both),
"same_dir_pear_qualified":True,"KEGG_tests_fixed":151,"KEGG_mapped_universe":4775,
"flavonoid":{r["Candidate_Set"]:{"hits":r["Genes_in_Pathway"],"FDR":r["BH_FDR_151_tested"]} for r in fig5},
"source_release":"2026-10-20", "warnings":["Biological independence of pooled primary apple sampling remains unresolved.","The v10 synteny rule is new and does not reproduce the original incomplete legacy scores.","Two PKS5-related loci drive the small-set flavonoid enrichment."]}
(O/"V10_FIGURE_TABLE_COUNTS.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2,ensure_ascii=False))
print("Supplementary rows:",{"S1":len(out),"S2":len(keep),"S3":len(ext),"S4":len(both),"S5":len(k2),"S6":len(p2),"S7":len(s7)})
