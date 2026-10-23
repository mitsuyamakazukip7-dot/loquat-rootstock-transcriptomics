#!/usr/bin/env python3
"""v10: fully specified cross-species local collinearity, independently of function.
Input: original strict 1:1:1 RBH triples, original GFF/feature_table, original expression table.
Primary definition: 10 neighboring strict RBH anchor genes in apple (5 each side);
for each target species consider homologs within 10 RBH-anchor ranks on the same
chromosome as the focal target; largest monotone (forward/reverse) subsequence
among valid homologs gives fraction of 10; final score is MIN(pear, loquat).
Threshold >=0.8 means at least eight locally collinear apple anchors in BOTH
pear and loquat. No annotations are consulted in filter/rank.
"""
import os
from pathlib import Path
from collections import defaultdict,Counter
import csv,gzip,re,math,json,bisect
ROOT=Path(os.environ.get("LOQUAT_PROJECT_ROOT",str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
OUT=ROOT/"11_methods_audit_20261020/v10_universal_synteny"
OUT.mkdir(parents=True,exist_ok=True)
def rows(p):
 with open(p,newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f))
original=rows(ROOT/"07_loquat_projection/cross_species/all_strict_ortholog_scored.csv")
legacy=rows(ROOT/"09_tables/Table_candidate_blind_ranking_all20795.csv")
legacy_lookup={r["apple_geneid"]:r for r in legacy}
coords={}
with gzip.open(ROOT/"03_reference/apple_GDDH13/GCF_002114115.1_ASM211411v1_feature_table.txt.gz","rt") as f:
 for r in csv.DictReader(f,delimiter="\t"):
  if r.get("# feature")=="gene" and r.get("GeneID"):
   coords["a",r["GeneID"]]=(r.get("genomic_accession") or r.get("chromosome"),int(r["start"]))
for sp,file in (("p","03_reference/pear_Cuiguan_GWHBAOS/GWHBAOS00000000.gff.gz"),("l","03_reference/loquat_TBY/GWHETRZ00000000.1.gff.gz")):
 with gzip.open(ROOT/file,"rt") as f:
  for l in f:
   x=l.rstrip("\n").split("\t")
   if len(x)!=9 or x[2]!="gene":continue
   m=re.search(r"(?:^|;)ID=([^;]+)",x[8])
   if m:coords[sp,m.group(1)]=(x[0],int(x[3]))
genespec={"a":"apple_geneid","p":"pear_evm","l":"loquat_eja"}
ordered={};index={}
for sp,key in genespec.items():
 d=defaultdict(list)
 for r in original:
  v=coords.get((sp,r[key]))
  if v:d[v[0]].append((v[1],r["apple_geneid"]))
 ordered[sp]={chrom:[i for _,i in sorted(entries)] for chrom,entries in d.items()}
 index[sp]={g:(chrom,idx) for chrom,vals in ordered[sp].items() for idx,g in enumerate(vals)}
def lis_len(v):
 tail=[]
 for x in v:
  p=bisect.bisect_left(tail,x)
  if p==len(tail):tail.append(x)
  else:tail[p]=x
 return len(tail)
def neighboring_apple(g,n=5):
 if g not in index["a"]:return []
 chrom,i=index["a"][g];x=ordered["a"][chrom]
 return x[max(0,i-n):i]+x[i+1:i+n+1]
def per_species(g,neighbors,sp,width=10):
 if g not in index[sp] or not neighbors:return 0,0,[]
 chrom,i=index[sp][g]
 # original array order is apple chromosome-order
 z=[index[sp][a][1] for a in neighbors if a in index[sp] and index[sp][a][0]==chrom and abs(index[sp][a][1]-i)<=width]
 if not z:return 0,0,[]
 m=max(lis_len(z),lis_len([-k for k in z]))
 return m/len(neighbors),m,z
def val(r,k):
 try: return float(r.get(k) or "nan")
 except: return float("nan")
def sig(r,k):return val(r,k)<0.05
def primary_class(r):
 if r["apple_core"]=="True":return "stage_consistent"
 if r["apple_stage_modulated"]=="True" and (sig(r,"apple_padj_s1") or sig(r,"apple_padj_s2")):return "interaction_associated"
 return "other"
def pear_ok(r):
 a,b=val(r,"apple_lfc_s1"),val(r,"apple_lfc_s2")
 if r["apple_core"]=="True": ref=(a+b)/2
 elif sig(r,"apple_padj_s1") and not sig(r,"apple_padj_s2"):ref=a
 elif sig(r,"apple_padj_s2") and not sig(r,"apple_padj_s1"):ref=b
 else:ref=(a+b)/2
 p=val(r,"pear_lfc")
 return math.isfinite(ref) and math.isfinite(p) and abs(p)>=1 and ref*p>0
result=[]
for r in original:
 g=r["apple_geneid"]
 ns=neighboring_apple(g)
 pear,ph,pz=per_species(g,ns,"p")
 loquat,lh,lz=per_species(g,ns,"l")
 cls=primary_class(r)
 is_eligible=cls!="other" and val(r,"loq_tpm_median")>=1 and pear_ok(r)
 ext=sig(r,"validation_new_padj") or sig(r,"validation_fruit_padj")
 both=str(r["validation_both_same_dir"]).lower()=="true"
 d={"apple_geneid":g,"pear_evm":r["pear_evm"],"loquat_eja":r["loquat_eja"],
    "class":cls,"n_reference_apple_anchors":len(ns),
    "pear_collinear_n":ph,"loquat_collinear_n":lh,
    "pear_collinear_fraction":round(pear,6),"loquat_collinear_fraction":round(loquat,6),
    "synteny_score_v10_min":round(min(pear,loquat),6),
    "synteny_score_v10_mean":round((pear+loquat)/2,6),
    "eligible_expression_direction":is_eligible,
    "external_any":ext,"external_both":both,
    "apple_padj_s1":r["apple_padj_s1"],"apple_padj_s2":r["apple_padj_s2"],
    "apple_interaction_padj":r["apple_interaction_padj"],
    "pear_lfc":r["pear_lfc"],"pear_qval":r["pear_qval"],
    "external_new_padj":r["validation_new_padj"],
    "external_fruit_padj":r["validation_fruit_padj"],
    "loquat_stem_median_TPM":r["loq_tpm_median"],
    "legacy_synteny_score":legacy_lookup.get(g,{}).get("synteny_score_blind","")}
 result.append(d)
names=list(result[0])
def save(name,rows):
 with open(OUT/name,"w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(rows)
save("all_20795.csv",result)
eligible=[r for r in result if r["eligible_expression_direction"]]
save("eligible_269.csv",eligible)
prior_ids={r["apple_geneid"] for r in rows(ROOT/"09_tables/submission_v9/TableS1_candidate_pool38.csv")}
summary={"algorithm":"At most 5 neighboring strict RBH anchors on each side in apple, target within 10 RBH ranks, monotone orientation-free LIS, min of pear and loquat fractions","n_strict_triples":len(result),"eligible_genes":len(eligible),"n_reference_apple_anchors_min":min(r["n_reference_apple_anchors"] for r in eligible)}
for th in (0.6,0.7,0.8,0.9,1.0):
 # Enforce the archived v10 minimum-coverage requirement, not just its fraction:
 # otherwise three boundary genes with only 1–3 anchors pass and inflate 190 to 193.
 pool=[r for r in eligible if r["n_reference_apple_anchors"]>=8 and r["synteny_score_v10_min"]>=th]
 ext=[r for r in pool if r["external_any"]]
 both=[r for r in pool if r["external_both"]]
 label=f"{th:.1f}"
 save(f"pool_{label}.csv",pool)
 save(f"external_{label}.csv",ext)
 save(f"both_{label}.csv",both)
 if th==0.8:
  save("s0.8_b0_pool.csv",pool)
  save("s0.8_b0_external.csv",ext)
  save("s0.8_b0_both.csv",both)
 summary[label]={"pool":len(pool),"external_any":len(ext),"external_both":len(both),
                 "intersect_legacy38":len({r["apple_geneid"] for r in pool}&prior_ids),
                 "legacy38_lost":len(prior_ids-{r["apple_geneid"] for r in pool}),
                 "flavonoid_annotated_candidates":sorted(set(r["apple_geneid"] for r in pool)&{"103437326","103421794","103443512"})}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2,ensure_ascii=False))
