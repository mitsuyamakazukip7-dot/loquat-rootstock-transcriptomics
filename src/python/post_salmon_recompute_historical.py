#!/usr/bin/env python3
"""Rebuild downstream pipeline in a new workspace; preserve original analysis untouched.
Runs only if all 51 independently downloaded/quantified samples validate.
"""
from pathlib import Path
import os, sys, csv, json, shutil, subprocess, time, datetime, hashlib
ORIG=Path(os.environ.get("LOQUAT_SOURCE_DIR", str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
BASE=Path(__file__).resolve().parent
WS=BASE/"recomputed_workspace"
LOG=BASE/"post_pipeline_log.txt"
def event(message):
 t=f"[{datetime.datetime.now(datetime.timezone.utc).isoformat()}] {message}"
 print(t,flush=True)
 with LOG.open("a") as f:f.write(t+"\n")
def cp(rel):
 src=ORIG/rel;dst=WS/rel
 if not src.is_file():raise FileNotFoundError(src)
 dst.parent.mkdir(exist_ok=True,parents=True)
 shutil.copy2(src,dst)
def patched(rel):
 src=ORIG/rel;dst=WS/rel
 dst.parent.mkdir(parents=True,exist_ok=True)
 s=src.read_text()
 s=s.replace("$HOME/Desktop/LOQUAT_DWARFING_PUBLIC",str(WS))
 s=s.replace("~/Desktop/LOQUAT_DWARFING_PUBLIC",str(WS))
 s=s.replace('"Desktop/LOQUAT_DWARFING_PUBLIC"',f'"{WS}"')
 s=s.replace("'Desktop/LOQUAT_DWARFING_PUBLIC'",f"'{WS}'")
 s=s.replace("Path.home()/'Desktop'/'LOQUAT_DWARFING_PUBLIC'",f"Path('{WS}')")
 # Ensure new pipeline reads a frozen 2026-10-20 mapping, not live KEGG
 if src.name=="reviewer_v8_pathway_label_permutation.R":
  start=s.index("links <- keggLink(")
  end=s.index("cc <- read.csv(",start)
  s=s[:start]+"""ttg <- read.csv('11_methods_audit_20261020/KEGG_frozen_TERM2GENE.csv',stringsAsFactors=FALSE)
sets <- split(as.character(ttg$gene),as.character(ttg$pathway))
sets <- lapply(sets,intersect,y=genes)
sets <- sets[lengths(sets)>=10 & lengths(sets)<=500]
"""+s[end:]
 if src.name=="apple_main36_deseq2.R":
  old="g<-gseKEGG(geneList=gl,organism='mdm',keyType='ncbi-geneid',minGSSize=10,maxGSSize=500,pvalueCutoff=1,verbose=FALSE)"
  new="kegg_ttg<-read.csv(file.path(root,'11_methods_audit_20261020','KEGG_frozen_TERM2GENE.csv')); g<-clusterProfiler::GSEA(geneList=gl,TERM2GENE=kegg_ttg[,c('pathway','gene')],minGSSize=10,maxGSSize=500,pvalueCutoff=1,verbose=FALSE,seed=TRUE)"
  if old not in s:raise AssertionError("Cannot patch live gseKEGG call")
  s=s.replace(old,new)
 if src.name=="v10_universal_synteny.py":
  s=s.replace('legacy_lookup[g].get("synteny_score_blind","")','legacy_lookup.get(g,{}).get("synteny_score_blind","")')
 if src.name=="run_pairwise_rbh.sh":
  # Avoid dependency on pandas: use Python standard-library CSV for pairwise RBH reduction.
  beg=s.index("import pandas as pd, os")
  end=s.index("\nPY",beg)
  impl="""import csv, os
root=os.path.expanduser('__RBH_ROOT__')
def read(path):
 with open(path,newline='') as f:
  r=csv.reader(f,delimiter='\t')
  return [dict(zip(('q','s','pid','len','qlen','slen','e','bits'),v)) for v in r]
for a,b in [('apple','pear'),('apple','loquat'),('pear','loquat')]:
 rev={r['q']:r['s'] for r in read(f'{root}/{b}_vs_{a}.tsv')}
 hits=[r for r in read(f'{root}/{a}_vs_{b}.tsv') if rev.get(r['s'])==r['q']]
 with open(f'{root}/{a}_{b}_RBH.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['q','s','pid','len','qlen','slen','e','bits'],delimiter='\t')
  w.writeheader();w.writerows(hits)
 print(a,b,len(hits))
""".replace('__RBH_ROOT__',str(WS/'05_ortholog/rbh'))
  s=s[:beg]+impl+s[end:]
 dst.write_text(s)
def run(cmd,stage):
 event("START "+stage)
 with (BASE/"post_commands.log").open("a") as l:
  l.write("\n\n### "+stage+": "+" ".join(map(str,cmd))+"\n");l.flush()
  p=subprocess.run(list(map(str,cmd)),cwd=WS,stdout=l,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError(f"{stage}: command exited {p.returncode}")
 event("PASS "+stage)
def check_sources():
 expected={"apple_main36":36,"apple2026":12,"loquat_TBY_stem":3}
 for name,n in expected.items():
  qs=list((BASE/"quant"/name).glob("*/quant.sf"))
  if len(qs)!=n:raise RuntimeError(f"{name} has {len(qs)} quant.sf but requires {n}")
def prepare():
 WS.mkdir(parents=True,exist_ok=True)
 ref=WS/"03_reference"
 if not ref.exists():ref.symlink_to(ORIG/"03_reference",target_is_directory=True)
 qa=WS/"04_quant"
 if not qa.exists():qa.symlink_to(BASE/"quant",target_is_directory=True)
 for sub in ["05_ortholog/rbh","06_meta_analysis/apple_main36","06_meta_analysis/apple2026","06_meta_analysis/pear2026/pathway","06_meta_analysis/reviewer_v8",
             "07_loquat_projection/cross_species","08_figures/apple_main36","08_figures/apple2026",
             "09_tables/submission_v9","11_methods_audit_20261020","12_submission_v10_assets"]:
  (WS/sub).mkdir(parents=True,exist_ok=True)
 inputs=[
 "01_metadata/apple_main36_samples.tsv","01_metadata/apple2026_samples.tsv","01_metadata/loquat_TBY_stem_samples.tsv",
 "06_meta_analysis/pear2026_deg_full.csv","06_meta_analysis/pear2026/pathway/pear_KEGG_enrichment.csv",
 "09_tables/Table_candidate_blind_ranking_all20795.csv",
 "09_tables/submission_v9/TableS1_candidate_pool38.csv","09_tables/submission_v9/TableS2_external_supported17.csv",
 "11_methods_audit_20261020/KEGG_frozen_TERM2GENE.csv",
 "11_methods_audit_20261020/KEGG_frozen_mapped_background.txt"]
 for x in inputs:cp(x)
 scripts=["scripts/apple_main36_deseq2.R","scripts/apple2026_deseq2.R","scripts/apple2026_concordance.R",
  "scripts/run_pairwise_rbh.sh","scripts/integrate_cross_species.py","scripts/reviewer_v8_pathway_label_permutation.R",
  "11_methods_audit_20261020/v10_universal_synteny.py","11_methods_audit_20261020/KEGG_v10_threshold_sensitivity.R",
  "11_methods_audit_20261020/KEGG_v10_PKS_leaveout.R","11_methods_audit_20261020/v10_package.py",
  "12_build_supp_v10.py","12_generate_4tables_v10.py","12_build_submission_tables_v10.py",
  "12_make_figs_v10.R"]
 for s in scripts:patched(s)
 # Recheck no unresolved static original path in scripts (the workspace path is permitted)
 for s in scripts:
  t=(WS/s).read_text()
  if str(ORIG) in t and str(WS) not in t:raise RuntimeError(f"Unpatched hardcoded source path {s}")
 event("Prepared isolated workspace")
def triples():
 base=WS/"05_ortholog/rbh"
 def load(filename):
  with (base/filename).open() as f:return list(csv.DictReader(f,delimiter="\t"))
 ap=load("apple_pear_RBH.tsv");al=load("apple_loquat_RBH.tsv")
 almap={x["q"]:x["s"] for x in al}
 out=base/"apple_pear_loquat_RBH_triples.tsv"
 with out.open("w",newline="") as f:
  w=csv.writer(f,delimiter="\t");w.writerow(["apple_protein","pear_protein","loquat_protein"])
  for r in ap:
   if r["q"] in almap:w.writerow([r["q"],r["s"],almap[r["q"]]])
 event(f"Built apple-centered protein RBH triples {sum(1 for _ in out.open())-1}; no pear-loquat RBH requirement")
def compare():
 summary={"datetime":datetime.datetime.now(datetime.timezone.utc).isoformat(),"frozen_KEGG_used":True}
 def src(rel):return WS/rel
 try:
  r=list(csv.DictReader(src("06_meta_analysis/apple_main36/DE_summary.csv").open()))
  summary["primary_DE"]={x[""]:int(x["FDR05"]) for x in r}
 except Exception as e:summary["primary_DE_error"]=str(e)
 try:
  r=list(csv.DictReader(src("06_meta_analysis/apple2026/DE_summary.csv").open()))
  summary["external_DE"]={x[""]:int(x["FDR05"]) for x in r}
 except Exception as e:summary["external_DE_error"]=str(e)
 try:
  r=json.loads(src("11_methods_audit_20261020/V10_RESULTS.json").read_text())
  summary["candidate_counts"]=r["counts"];summary["flavonoid_FDR"]=r["flavonoid_fixed_151_pathway_FDR"]
 except Exception as e:summary["candidate_counts_error"]=str(e)
 (BASE/"FULL_RECOMPUTATION_RESULTS.json").write_text(json.dumps(summary,indent=2)+"\n")
 event("Wrote result comparison "+str(summary))
def main():
 check_sources();prepare()
 steps=[
 (["Rscript","scripts/apple_main36_deseq2.R"],"main_de_and_frozen_kegg"),
 (["Rscript","scripts/apple2026_deseq2.R"],"external_de"),
 (["Rscript","scripts/apple2026_concordance.R"],"external_concordance"),
 (["zsh","scripts/run_pairwise_rbh.sh"],"diamond_pairwise_RBH")]
 try:
  for cmd,stage in steps:run(cmd,stage)
  triples()
  run(["python3","scripts/integrate_cross_species.py"],"integrate_three_species")
  run(["Rscript","scripts/reviewer_v8_pathway_label_permutation.R"],"permutation_fixed_KEGG")
  run(["python3","11_methods_audit_20261020/v10_universal_synteny.py"],"collinearity")
  run(["Rscript","11_methods_audit_20261020/KEGG_v10_threshold_sensitivity.R"],"KEGG_frozen_sensitivity")
  run(["Rscript","11_methods_audit_20261020/KEGG_v10_PKS_leaveout.R"],"PKS_leaveout")
  run(["python3","11_methods_audit_20261020/v10_package.py"],"v10_QA_package")
  run(["python3","12_build_supp_v10.py"],"supplementary_tables")
  compare()
  event("DOWNSTREAM RECOMPUTATION FINISHED")
 except Exception as exc:
  compare()
  event("FAILED "+str(exc));sys.exit(1)
if __name__=="__main__":main()
