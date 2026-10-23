#!/usr/bin/env python3
"""Independent public FASTQ -> Salmon quant.sf runner. Never edits original analysis."""
import os
from pathlib import Path
import argparse, csv, json, datetime, subprocess, shutil, sys, os, time
ORIG=Path(os.environ.get("LOQUAT_SOURCE_DIR", str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
BASE=Path(__file__).resolve().parent
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def say(s):print(f"[{now()}] {s}",flush=True)
def rows(file,sep="\t"):
 with open(file,newline="") as f:return list(csv.DictReader(f,delimiter=sep))
def samples():
 s=[]
 for r in rows(ORIG/"01_metadata/apple_main36_samples.tsv"):s.append(dict(dataset="apple_main36",name=r["sample"],run=r["Run"],ref="apple_GDDH13"))
 for r in rows(ORIG/"01_metadata/apple2026_samples.tsv"):s.append(dict(dataset="apple2026",name=r["sample"],run=r["run"],ref="apple_GDDH13"))
 for r in rows(ORIG/"01_metadata/loquat_TBY_stem_samples.tsv"):s.append(dict(dataset="loquat_TBY_stem",name=r["sample"],run=r["run"],ref="loquat_TBY"))
 assert len(s)==51 and len({x["run"] for x in s})==51
 return s
def execute(cmd,log):
 cmd=list(map(str,cmd))
 say("RUN: "+" ".join(cmd))
 with log.open("a") as f:
  f.write("\n["+now()+"] "+" ".join(cmd)+"\n");f.flush()
  p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError(f"{cmd[0]} exit {p.returncode}; details: {log}")
def quant_ok(out):
 q=out/"quant.sf";m=out/"aux_info/meta_info.json"
 if not (q.is_file() and m.is_file()):return False
 try:
  meta=json.loads(m.read_text())
  if meta.get("num_mapped",0)<100:return False
  total=0.;n=0
  for row in rows(q):
   total+=float(row["TPM"]);n+=1
  return n>10000 and 999990<total<1000010
 except: return False
def process(s):
 dataset,name,run=s["dataset"],s["name"],s["run"]
 log=BASE/"logs"/dataset/(name+".log");log.parent.mkdir(parents=True,exist_ok=True)
 qdir=BASE/"quant"/dataset/name
 status=BASE/"status"/dataset/(name+".json");status.parent.mkdir(parents=True,exist_ok=True)
 obj=dict(dataset=dataset,sample=name,run=run,started=now(),status="RUNNING")
 def save():status.write_text(json.dumps(obj,indent=2)+"\n")
 save()
 if quant_ok(qdir):
  obj.update(status="SKIP_VERIFIED",finished=now());save();say(f"SKIP_VERIFIED {name}");return
 if qdir.exists():raise RuntimeError(f"Existing invalid/incomplete quant requires inspection: {qdir}")
 free=shutil.disk_usage(BASE).free/1024**3
 if free<45:raise RuntimeError(f"disk reserve violated: only {free:.1f} GiB free")
 scratch=BASE/"scratch"/run;scratch.mkdir(parents=True,exist_ok=True)
 try:
  # Prefetch archive to isolated scratch, validate, create temporary paired FASTQ
  execute(["prefetch","--max-size","30G","-O",scratch/"sra",run],log)
  sr=list((scratch/"sra").rglob("*.sra"))
  if len(sr)!=1:raise RuntimeError(f"Expected 1 SRA file, found {sr}")
  execute(["vdb-validate",sr[0]],log)
  fastq=scratch/"fastq";fastq.mkdir(exist_ok=True)
  tmp=scratch/"tmp";tmp.mkdir(exist_ok=True)
  execute(["fasterq-dump","--split-files","--threads","2","--temp",tmp,"--outdir",fastq,sr[0]],log)
  a,b=fastq/(run+"_1.fastq"),fastq/(run+"_2.fastq")
  if not a.is_file() or not b.is_file() or min(a.stat().st_size,b.stat().st_size)<1024:raise RuntimeError("missing paired FASTQ")
  obj.update(sra_bytes=sr[0].stat().st_size,fastq_bytes=[a.stat().st_size,b.stat().st_size]);save()
  idx=ORIG/"03_reference"/s["ref"]/"salmon_index"
  qdir.parent.mkdir(parents=True,exist_ok=True)
  execute(["salmon","quant","-i",idx,"-l","A","-1",a,"-2",b,"-p","4","--validateMappings","--gcBias","-o",qdir],log)
  if not quant_ok(qdir):raise RuntimeError(f"Failed quant.sf structure/TPM validation: {qdir}")
  meta=json.loads((qdir/"aux_info/meta_info.json").read_text())
  obj.update(status="PASS",mapped_percent=meta.get("percent_mapped"),processed=meta.get("num_processed"),finished=now());save()
  shutil.rmtree(scratch)
  say(f"PASS {dataset}/{name} mapped={obj['mapped_percent']}")
 except Exception as e:
  obj.update(status="FAILED",error=str(e),failed_at=now());save()
  raise
def main():
 parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group(required=True);g.add_argument("--pilot",action="store_true");g.add_argument("--all",action="store_true");args=parser.parse_args()
 BASE.mkdir(parents=True,exist_ok=True)
 seq=samples()[:1] if args.pilot else samples()
 progress=dict(mode="pilot" if args.pilot else "all",total=len(seq),completed=0,started=now())
 for i,s in enumerate(seq,1):
  say(f"START {i}/{len(seq)} {s['dataset']}/{s['name']} {s['run']}")
  try:process(s)
  except Exception as e:
   progress.update(error=str(e),failed=s,stopped=now())
   (BASE/"progress.json").write_text(json.dumps(progress,indent=2)+"\n")
   say(f"STOP: {e}");sys.exit(1)
  progress.update(completed=i,latest=s,updated=now())
  (BASE/"progress.json").write_text(json.dumps(progress,indent=2)+"\n")
 say("ALL REQUESTED FASTQ QUANTIFICATIONS COMPLETE")
if __name__=="__main__":main()
