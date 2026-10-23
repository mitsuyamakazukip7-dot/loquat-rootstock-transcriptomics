#!/usr/bin/env python3
"""Portable independent post-run checks (do not use checks to generate results)."""
from pathlib import Path
import csv,json,math
ROOT=Path(__file__).resolve().parent
def rows(name):
 with (ROOT/name).open(newline="") as f:return list(csv.DictReader(f))
x=rows("PERMUTATION_REPRODUCTION_COMPARISON.csv")
lookup={r["subset"]:r for r in x}
for key, n,obs,tail,p in [
 ("historical_selected37",37,27,158,.395),
 ("frozen_full_selected38",38,28,158,.395),
 ("all147",147,86,168,.420)]:
 r=lookup[key]
 assert int(r["n_pathways"])==n,(key,"n")
 assert int(r["observed_reversals"])==obs,(key,"obs")
 assert int(r["tail_ge"])==tail,(key,"tail")
 assert math.isclose(float(r["p_ge"]),p,abs_tol=1e-12),(key,"p")
q=rows("INDEPENDENT_VS_COUPLED_PERMUTATION.csv")
assert len(q)==3
m=rows("KEGG_MDM_FULL_TERM2GENE_20261022.csv")
assert len(m)==21003
assert len({r["pathway"] for r in m})==163
a=rows("GSEA_FULL_FROZEN_STAGE_COMPARISON.csv")
assert len(a)==147
assert sum(r["sharedFDR10"].lower()=="true" for r in a)==38
print("PASS portable frozen GSEA/permutation: 37/27/.395; 38/28/.395; 147/86/.420")
(ROOT/"PORTABLE_RECHECK_PASS.json").write_text(json.dumps({"passed":True,"historical37":{"observed":27,"p":0.395},"deterministic38":{"observed":28,"p":.395},"all147":{"observed":86,"p":.42}},indent=2))
