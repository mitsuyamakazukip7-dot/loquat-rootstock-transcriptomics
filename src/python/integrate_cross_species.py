import os
from pathlib import Path
import csv, gzip, re, math, statistics
from collections import Counter, defaultdict

ROOT = Path(os.environ.get("LOQUAT_PROJECT_ROOT",str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
OUT = ROOT/"07_loquat_projection"/"cross_species"
OUT.mkdir(parents=True, exist_ok=True)

def num(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

def sig(x, cutoff=0.05):
    v=num(x); return math.isfinite(v) and v < cutoff

def rows(path, sep=','):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter=sep))

def index(rs, key):
    return {str(r.get(key,'')): r for r in rs if str(r.get(key,''))}

def sign(x):
    x=num(x)
    return 1 if x>0 else (-1 if x<0 else 0)
# Apple protein accession -> GeneID
apple_p2g = {}
ft = ROOT/"03_reference/apple_GDDH13/GCF_002114115.1_ASM211411v1_feature_table.txt.gz"
with gzip.open(ft, 'rt', encoding='utf-8') as f:
    r = csv.DictReader(f, delimiter='\t')
    for x in r:
        p=x.get('product_accession',''); g=x.get('GeneID','')
        if p and g: apple_p2g[p]=g

def parse_headers(path, species):
    p2g, tx2g = {}, {}
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if not line.startswith('>'): continue
            pid=line[1:].split()[0]
            fields=dict(re.findall(r'(\w+)=([^\t ]+)', line))
            if species=='pear': p2g[pid]=fields.get('OriGeneID','')
            else:
                g=fields.get('OriGeneID',''); tx=fields.get('mRNA','')
                p2g[pid]=g
                if tx and g: tx2g[tx]=g
    return p2g, tx2g

pear_p2g,_ = parse_headers(ROOT/"03_reference/pear_Cuiguan_GWHBAOS/GWHBAOS00000000.Protein.faa.gz",'pear')
loq_p2g,loq_tx2g = parse_headers(ROOT/"03_reference/loquat_TBY/GWHETRZ00000000.1.Protein.faa.gz",'loquat')
# Protein RBH triples -> strict gene-level 1:1:1 triples
tri = rows(ROOT/"05_ortholog/rbh/apple_pear_loquat_RBH_triples.tsv", '\t')
gene_triples = {}
for t in tri:
    a=apple_p2g.get(t['apple_protein'],'')
    p=pear_p2g.get(t['pear_protein'],'')
    l=loq_p2g.get(t['loquat_protein'],'')
    if not (a and p and l): continue
    gene_triples[(a,p,l)] = (t['apple_protein'],t['pear_protein'],t['loquat_protein'])
ac=Counter(k[0] for k in gene_triples); pc=Counter(k[1] for k in gene_triples); lc=Counter(k[2] for k in gene_triples)
strict=[]
for (a,p,l),(ap,pp,lp) in gene_triples.items():
    if ac[a]==pc[p]==lc[l]==1:
        strict.append({'apple_geneid':a,'pear_evm':p,'loquat_eja':l,
                       'apple_protein':ap,'pear_protein':pp,'loquat_protein':lp})

# Loquat stem transcript TPM -> Eja gene TPM
samples=['Stem-1','Stem-2','Stem-3']; gtpm=defaultdict(lambda:[0.0,0.0,0.0])
for i,s in enumerate(samples):
    q=ROOT/f"04_quant/loquat_TBY_stem/{s}/quant.sf"
    for r in rows(q,'\t'):
        g=loq_tx2g.get(r['Name'],'')
        if g: gtpm[g][i]+=num(r['TPM'])
loq_expr={}
for g,v in gtpm.items():
    m=sum(v)/3; med=statistics.median(v); mn=min(v)
    sd=statistics.stdev(v) if len(v)>1 else 0.0
    loq_expr[g]={'loq_tpm1':v[0],'loq_tpm2':v[1],'loq_tpm3':v[2],
                 'loq_tpm_mean':m,'loq_tpm_median':med,'loq_tpm_min':mn,
                 'loq_tpm_cv':(sd/m if m>0 else float('nan'))}

# Analysis tables
base=ROOT/"06_meta_analysis"
a1=index(rows(base/"apple_main36/DE_stage1_M9_vs_MM106.csv"),'GENEID')
a2=index(rows(base/"apple_main36/DE_stage2_M9_vs_MM106.csv"),'GENEID')
ai=index(rows(base/"apple_main36/DE_stage_interaction.csv"),'GENEID')
v1=index(rows(base/"apple2026/DE_new_shoot_PYGD076_vs_M9T337.csv"),'GENEID')
v2=index(rows(base/"apple2026/DE_fruiting_spur_PYGD076_vs_M9T337.csv"),'GENEID')
pear=index(rows(base/"pear2026_deg_full.csv"),'gene_id')

patterns={
 'gibberellin':r'gibberellin|GA2-oxidase|GA20-oxidase|GA3-oxidase',
 'auxin':r'auxin|\bARF\b|\bPIN\b|\bIAA\b|\bSAUR\b|\bGH3\b',
 'brassinosteroid':r'brassinosteroid|\bBRI1\b|\bBZR|\bBES1\b|\bDWF4\b',
 'cell_wall':r'cell wall|cellulose|pectin|expansin|xyloglucan|lignin|polygalacturonase|pectate',
 'vascular':r'xylem|phloem|cambium|vascular|secondary wall|\bVND\b|\bNST\b|\bSND\b',
 'sugar_T6P':r'trehalose|T6P|SnRK1|sucrose|starch|trehalose-phosphate synthase|trehalose phosphatase'
}
allrows=[]
for t in strict:
    g=t['apple_geneid']; p=t['pear_evm']; l=t['loquat_eja']
    r={**t}; A1=a1.get(g,{}); A2=a2.get(g,{}); INT=ai.get(g,{})
    P=pear.get(p,{}); V1=v1.get(g,{}); V2=v2.get(g,{})
    l1=num(A1.get('log2FoldChange')); l2=num(A2.get('log2FoldChange'))
    s1=sig(A1.get('padj')); s2=sig(A2.get('padj')); si=sig(INT.get('padj'))
    core=s1 and s2 and sign(l1)==sign(l2) and sign(l1)!=0
    anchor=(l1+l2)/2 if core else (l1 if s1 and not s2 else (l2 if s2 and not s1 else (l1+l2)/2))
    plfc=num(P.get('log2(fc)')); ps=sig(P.get('qval')) and math.isfinite(plfc) and abs(plfc)>=1
    pc=ps and sign(anchor)==sign(plfc) and sign(anchor)!=0
    val1=sig(V1.get('padj')); val2=sig(V2.get('padj'))
    valboth=val1 and val2 and sign(V1.get('log2FoldChange'))==sign(V2.get('log2FoldChange'))
    e=loq_expr.get(l,{'loq_tpm_mean':0,'loq_tpm_median':0,'loq_tpm_min':0,'loq_tpm_cv':float('nan'),
                      'loq_tpm1':0,'loq_tpm2':0,'loq_tpm3':0})
    desc=' '.join([A1.get('symbol',''),A1.get('name',''),P.get('Description',''),P.get('KEGG','')])
    cats=[k for k,pat in patterns.items() if re.search(pat,desc,re.I)]
    score=(4 if core else (2 if (s1 or s2) else 0)) + (1 if si else 0) + (2 if ps else 0) + (2 if pc else 0)
    score += (1 if (val1 or val2) else 0) + (1 if valboth else 0) + (1 if e['loq_tpm_median']>=1 else 0)
    score += (1 if e['loq_tpm_median']>=5 else 0) + (1 if e['loq_tpm_mean']>=1 and e['loq_tpm_cv']<=0.5 else 0) + (1 if cats else 0)
    r.update({'apple_symbol':A1.get('symbol',''),'apple_name':A1.get('name',''),
              'apple_lfc_s1':l1,'apple_padj_s1':A1.get('padj',''),'apple_lfc_s2':l2,'apple_padj_s2':A2.get('padj',''),
              'apple_interaction_padj':INT.get('padj',''),'apple_core':core,'apple_stage_modulated':si,
              'pear_lfc':plfc,'pear_qval':P.get('qval',''),'pear_regulation':P.get('regulation',''),
              'pear_description':P.get('Description',''),'pear_kegg':P.get('KEGG',''),'pear_sig':ps,'pear_concordant':pc,
              'validation_new_padj':V1.get('padj',''),'validation_fruit_padj':V2.get('padj',''),
              'validation_any':val1 or val2,'validation_both_same_dir':valboth,
              'categories':'|'.join(cats),'score':score,**e})
    if core and ps and pc and e['loq_tpm_median']>=1:
        r['class']='conserved_core'
    elif (s1 or s2) and si and ps and pc and e['loq_tpm_median']>=1:
        r['class']='stage_modulated_conserved'
    elif (s1 or s2) and e['loq_tpm_median']>=1:
        r['class']='apple_rootstock_loquat_expressed'
    else: r['class']='other'
    allrows.append(r)

allrows.sort(key=lambda x:(-x['score'],-num(x['loq_tpm_median'])))
fields=list(allrows[0].keys())
with open(OUT/'all_strict_ortholog_scored.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(allrows)
def write_subset(name, pred):
    sub=[r for r in allrows if pred(r)]
    with open(OUT/name,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(sub)
    return len(sub)

ncore=write_subset('conserved_core_candidates.csv',lambda r:r['class']=='conserved_core')
nstage=write_subset('stage_modulated_conserved_candidates.csv',lambda r:r['class']=='stage_modulated_conserved')
ntop=write_subset('top_candidates_score_ge8.csv',lambda r:r['score']>=8 and num(r['loq_tpm_median'])>=1)

# Pathway-level Apple-Pear concordance
ak=rows(base/"apple_main36/cross_stage_KEGG_concordance.csv")
pk=rows(base/"pear2026/pathway/pear_KEGG_enrichment.csv")
pkidx={re.search(r'(\d{5})$',r.get('code','') or r.get('Pathway_ID','')).group(1):r
       for r in pk if re.search(r'(\d{5})$',r.get('code','') or r.get('Pathway_ID',''))}
pout=[]
for a in ak:
    m=re.search(r'(\d{5})$',a['ID'])
    if not m: continue
    p=pkidx.get(m.group(1),{})
    aok=(a.get('same_direction')=='TRUE' and num(a.get('p.adjust_s1'))<.1 and num(a.get('p.adjust_s2'))<.1)
    pn=num(p.get('ZScore')); an=(num(a.get('NES_s1'))+num(a.get('NES_s2')))/2
    pq=num(p.get('Q.value')); psig=math.isfinite(pq) and pq<.05
    pout.append({'code':m.group(1),'apple_description':a.get('Description',''),'apple_NES_mean':an,
                 'apple_shared_FDR10':aok,'pear_pathway':p.get('Pathway_Name',''),'pear_ZScore':pn,
                 'pear_Q':p.get('Q.value',''),'pear_sig':psig,
                 'direction_concordant':aok and psig and sign(an)==sign(pn) and sign(an)!=0})
with open(OUT/'cross_species_pathway_concordance.csv','w',newline='') as f:
    pf=list(pout[0].keys()); w=csv.DictWriter(f,fieldnames=pf); w.writeheader(); w.writerows(pout)

summary={
 'protein_RBH_triples':len(tri),'distinct_gene_triples':len(gene_triples),'strict_gene_1to1to1':len(strict),
 'conserved_core_candidates':ncore,'stage_modulated_conserved_candidates':nstage,'score_ge8_loquat_expr':ntop,
 'apple_shared_KEGG_FDR10':sum(r['apple_shared_FDR10'] for r in pout),
 'apple_pear_pathway_sig_concordant':sum(r['direction_concordant'] for r in pout)
}
with open(OUT/'integration_summary.tsv','w') as f:
    for k,v in summary.items(): f.write(f'{k}\t{v}\n')

print('SUMMARY')
for k,v in summary.items(): print(k,v)
print('TOP15')
for r in allrows[:15]:
    print(r['score'],r['class'],r['apple_geneid'],r['apple_name'],r['pear_evm'],
          round(num(r['pear_lfc']),3),round(num(r['loq_tpm_median']),3),r['loquat_eja'],r['categories'])
