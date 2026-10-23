#!/bin/zsh
set -euo pipefail
ROOT="${LOQUAT_PROJECT_ROOT:-$HOME/Desktop/LOQUAT_DWARFING_PUBLIC}"
REF="$ROOT/03_reference"
OUT="$ROOT/05_ortholog/rbh"
LOG="$ROOT/logs/rbh"
mkdir -p "$OUT" "$LOG"
APPLE="$OUT/apple.faa"; PEAR="$OUT/pear.faa"; LOQUAT="$OUT/loquat.faa"
python3 - <<'PY'
import gzip,shutil,os
root=os.path.expanduser(os.environ.get('LOQUAT_PROJECT_ROOT','~/Desktop/LOQUAT_DWARFING_PUBLIC'))
items=[('03_reference/apple_GDDH13/GCF_002114115.1_ASM211411v1_protein.faa.gz','05_ortholog/rbh/apple.faa'),('03_reference/pear_Cuiguan_GWHBAOS/GWHBAOS00000000.Protein.faa.gz','05_ortholog/rbh/pear.faa'),('03_reference/loquat_TBY/GWHETRZ00000000.1.Protein.faa.gz','05_ortholog/rbh/loquat.faa')]
for src,dst in items:
 s=os.path.join(root,src); d=os.path.join(root,dst)
 if not os.path.exists(d):
  with gzip.open(s,'rb') as i, open(d,'wb') as o: shutil.copyfileobj(i,o)
PY
for sp in apple pear loquat; do
  [ -s "$OUT/${sp}.dmnd" ] || diamond makedb --in "$OUT/${sp}.faa" -d "$OUT/${sp}" > "$LOG/${sp}_makedb.log" 2>&1
done
pairs=("apple pear" "pear apple" "apple loquat" "loquat apple" "pear loquat" "loquat pear")
for p in $pairs; do
  q=${p%% *}; s=${p##* }
  out="$OUT/${q}_vs_${s}.tsv"
  [ -s "$out" ] && continue
  diamond blastp -q "$OUT/${q}.faa" -d "$OUT/${s}" -o "$out" \
    -f 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --sensitive --max-target-seqs 1 --evalue 1e-5 --query-cover 50 --subject-cover 50 -p 2 \
    > "$LOG/${q}_vs_${s}.log" 2>&1
done
python3 - <<'PY'
import pandas as pd, os
root=os.path.join(os.path.expanduser(os.environ.get('LOQUAT_PROJECT_ROOT','~/Desktop/LOQUAT_DWARFING_PUBLIC')),'05_ortholog/rbh')
for a,b in [('apple','pear'),('apple','loquat'),('pear','loquat')]:
 f1=pd.read_csv(f'{root}/{a}_vs_{b}.tsv',sep='\t',header=None,names=['q','s','pid','len','qlen','slen','e','bits'])
 f2=pd.read_csv(f'{root}/{b}_vs_{a}.tsv',sep='\t',header=None,names=['q','s','pid','len','qlen','slen','e','bits'])
 rev=dict(zip(f2.q,f2.s)); r=f1[f1.apply(lambda x: rev.get(x.s)==x.q,axis=1)]
 r.to_csv(f'{root}/{a}_{b}_RBH.tsv',sep='\t',index=False)
 print(a,b,len(r))
PY
