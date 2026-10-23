#!/usr/bin/env python3
"""Recreate the current v10 Fig1 and SuppFig1 on Mac from the same v10 data.
The rendering logic mirrors the previously approved three-panel / six-locus figures.
"""
import numpy as np,json,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Patch
from matplotlib.lines import Line2D
from pathlib import Path
import os
P=Path(os.environ.get("LOQUAT_PROJECT_ROOT",str(Path.home()/"Desktop/LOQUAT_DWARFING_PUBLIC"))).expanduser().resolve()
OUT=Path(os.environ.get("LOQUAT_FIGURE_OUT",str(P/"12_submission_v10_assets"/"LOCAL_FIG1_S1_RENDER"))).expanduser().resolve()
OUT.mkdir(exist_ok=True,parents=True)
C={"ink":"#26363B","sub":"#68767B","teal":"#8BD1C8","mint":"#B3DDD2",
"teal_deep":"#4C9E96","pink":"#F2BCD4","gold":"#D8BC72","cream":"#FFF4B3",
"pale":"#F4F8F7","stroke":"#BACBC9","gray":"#C7CED1","brown":"#9A765C","focal":"#D66F5B"}
plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
def box(ax,x,y,w,h,fc,ec=None,r=.022,lw=1.0):
 bb=FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.008,rounding_size={r}",
       facecolor=fc,edgecolor=ec or C["stroke"],linewidth=lw)
 ax.add_patch(bb);return bb
def txt(ax,x,y,s,size=11,color=None,weight="normal",ha="center",va="center",**kw):
 return ax.text(x,y,s,fontsize=size,color=color or C["ink"],
                fontweight=weight,ha=ha,va=va,**kw)
def arr(ax,xy1,xy2,ca="#93ADA9",ms=11,lw=1.4):
 ax.add_patch(FancyArrowPatch(xy1,xy2,arrowstyle="-|>",mutation_scale=ms,color=ca,linewidth=lw))
def savefig(fig,stem):
 for ext in ("png","pdf","svg"):
  fig.savefig(OUT/(stem+"."+ext),dpi=300 if ext=="png" else 150,
      facecolor="white",bbox_inches="tight",pad_inches=.15)

# ---------------- FIGURE 1 -----------------
fig=plt.figure(figsize=(15.2,12.5),facecolor="white")
a=fig.add_axes([.055,.735,.89,.245]);b=fig.add_axes([.055,.505,.89,.20])
c=fig.add_axes([.055,.035,.89,.445])
for ax in (a,b,c):
 ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis("off")
txt(a,0,1,"a",size=16,weight="bold",ha="left",va="top")
txt(a,.037,1,"Rootstock-associated differences in tree architecture",size=15,weight="bold",ha="left",va="top")
def tree(ax,cx,scale,name,large):
 base=.21
 ax.add_line(Line2D([cx,cx],[base,base+.34*scale],lw=9.5*scale,color=C["brown"],solid_capstyle="round"))
 for dx in [-.11,.11]:
  ax.add_line(Line2D([cx,cx+dx*scale],[base+.18*scale,base+.38*scale],lw=4*scale,color=C["brown"],solid_capstyle="round"))
 for dx,dy,r in [(-.13,.42,.09),(0,.51,.095),(.13,.42,.09),(-.055,.35,.095),(.065,.35,.09)]:
  ax.add_patch(Circle((cx+dx*scale,base+dy*scale),r*scale,
     facecolor=C["mint"] if large else C["teal"],edgecolor="white",lw=.7))
 ax.add_line(Line2D([cx,cx],[base-.045,base-.095],lw=8,color=C["brown"],solid_capstyle="butt"))
 ax.add_patch(Circle((cx,base-.033),.029,facecolor="white",edgecolor=C["gold"],lw=2))
 txt(ax,cx,.04,name,size=10.5,weight="bold")
tree(a,.23,.92,"Vigorous / standard",True)
tree(a,.74,.62,"Dwarfing rootstock",False)
box(a,.43,.27,.15,.34,"#FBFBFA",r=.025)
txt(a,.505,.5,"Same scion",size=11,weight="bold")
txt(a,.505,.405,"Different\nrootstocks",size=10)
txt(a,.505,.32,"Growth response",size=9,color=C["sub"])
arr(a,(.59,.45),(.63,.45),C["gold"],13)
txt(b,0,1,"b",size=16,weight="bold",ha="left",va="top")
txt(b,.037,1,"Orchard-level motivation",size=15,weight="bold",ha="left",va="top")
for i,(label,col) in enumerate(zip(
 ["Compact canopy","Higher planting density","Easier pruning / harvest","Simpler canopy management"],
 [C["mint"],C["teal"],C["pink"],C["cream"]])):
 x=.015+i*.25
 box(b,x,.12,.225,.62,col,r=.025,ec="white")
 txt(b,x+.1125,.42,label,size=10.5,weight="bold")
txt(c,0,1,"c",size=16,weight="bold",ha="left",va="top")
txt(c,.037,1,"Cross-species evidence framework",size=15,weight="bold",ha="left",va="top")
inputs=[(.01,.83,.245,.105,C["pale"],"Apple\nM9 vs MM106; two sampling contexts"),
        (.285,.83,.245,.105,"#FFF9DC","Pear\nDwarf germplasm vs Cuiguan"),
        (.56,.83,.245,.105,"#FFF0F5","Loquat\nStem expression + genome")]
for x,y,w,h,fc,lab in inputs:
 box(c,x,y,w,h,fc,r=.019);txt(c,x+w/2,y+h/2,lab,size=9.7,weight="medium")
box(c,.235,.66,.39,.10,"#EDF4F2",r=.025)
txt(c,.43,.71,"20,795 apple-centered RBH correspondences",size=13,weight="bold")
for x,y,w,h,fc,lab in inputs:arr(c,(x+w/2,y),(.43,.76),lw=.9,ms=7)
box(c,.24,.49,.38,.105,C["mint"],r=.02)
txt(c,.43,.553,"269 expression-supported triples",size=11.7,weight="bold")
txt(c,.43,.515,"23 stage-consistent  +  246 interaction-associated",size=9,color=C["sub"])
arr(c,(.43,.66),(.43,.6))
box(c,.24,.33,.38,.095,C["teal"],r=.02)
txt(c,.43,.385,"190 loci with conserved local gene order",size=11,weight="bold")
txt(c,.43,.352,"Min. pear/loquat score >= 0.80",size=9)
arr(c,(.43,.49),(.43,.43))
box(c,.24,.175,.38,.088,C["cream"],r=.02)
txt(c,.43,.218,"24 with external apple expression support",size=10.9,weight="bold")
arr(c,(.43,.33),(.43,.267))
box(c,.67,.195,.3,.115,"#FFF9DC",r=.025)
txt(c,.82,.267,"External apple dataset",size=10.8,weight="bold")
txt(c,.82,.229,"New-shoot leaves / fruiting-spur leaves",size=9,color=C["sub"])
arr(c,(.67,.25),(.626,.22),ca=C["gold"],ms=9)
box(c,.24,.04,.38,.08,"#E3F1EF",r=.02)
txt(c,.43,.080,"5 supported in both leaf contexts",size=11,weight="bold")
arr(c,(.43,.175),(.43,.124))
txt(c,.82,.1,"Functional annotation only after selection",size=9.1,color=C["sub"])
savefig(fig,"Fig1_V10_Overview");plt.close(fig)

# ---------------- SUPPLEMENTARY FIGURE 1 ------------------
obj=json.loads((P/"12_submission_v10_assets/SuppFig1_V10_6_loci_anchors.json").read_text())
assert len(obj["loci"])==6
fig,axs=plt.subplots(3,2,figsize=(15.2,13.1),sharex=True,sharey=True)
fig.patch.set_facecolor("white")
for i,(ax,rec) in enumerate(zip(axs.flat,obj["loci"])):
 ax.set_facecolor("white");ax.set_xlim(-6.7,6.0);ax.set_ylim(-.55,2.57)
 ax.spines[:].set_visible(False)
 ax.set_xticks(range(-5,6));ax.tick_params(axis="x",length=0,labelsize=8,pad=4,colors=C["sub"])
 ax.tick_params(axis="y",left=False,labelleft=False)
 tracks=rec["tracks"];ypos={"Apple":2,"Pear":1,"Loquat":0}
 focal=rec["apple_gene"]
 anchors={z["anchor"] for z in tracks[0]["genes"] if z["offset"]!=0}
 raw={t["species"]:{g["anchor"]:g["offset"] for g in t["genes"]} for t in tracks}
 polarity={"Apple":1}
 for target in ["Pear","Loquat"]:
  keys=[g for g in raw["Apple"] if g!=focal and g in raw[target]]
  xx=np.asarray([raw["Apple"][g] for g in keys],dtype=float)
  yy=np.asarray([raw[target][g] for g in keys],dtype=float)
  polarity[target]=(-1 if len(keys)>=2 and np.corrcoef(xx,yy)[0,1]<0 else 1)
 mapped={sp:{g:v*polarity[sp] for g,v in d.items()} for sp,d in raw.items()}
 for t in tracks:
  sp=t["species"];y=ypos[sp]
  ax.hlines(y,-5.52,5.52,color="#EFF1F0",lw=2.1,zorder=0)
  ax.text(-5.8,y,sp+(" *" if polarity.get(sp,1)<0 else ""),ha="right",va="center",fontsize=9.4,color=C["ink"])
 for s1,s2 in [("Apple","Pear"),("Pear","Loquat")]:
  m1,m2=mapped[s1],mapped[s2]
  for g in set(m1)&set(m2):
   col=C["focal"] if g==focal else C["teal_deep"] if g in anchors else C["stroke"]
   alpha=.72 if g==focal else .28 if g in anchors else .35
   wid=1.5 if g==focal else .7 if g in anchors else .5
   ax.plot([m1[g],m2[g]],[ypos[s1],ypos[s2]],color=col,alpha=alpha,lw=wid,zorder=1)
 for t in tracks:
  y=ypos[t["species"]]
  for g in t["genes"]:
   foc=g["focal"]
   col=C["focal"] if foc else C["mint"] if g["anchor"] in anchors else C["gray"]
   x=mapped[t["species"]][g["anchor"]]
   ax.add_patch(Rectangle((x-.31,y-.102),.62,.205,facecolor=col,
                  edgecolor="white" if not foc else "#A45145",lw=1 if foc else .35,zorder=3))
 ax.set_title(f"{chr(97+i)}  {rec['name']}  |  {rec['loquat_gene']}",
              fontweight="bold",loc="left",fontsize=12.2,color=C["ink"],pad=17)
 ax.text(0,2.37,f"Pear {float(rec['pear_score']):.2f}    Loquat {float(rec['loquat_score']):.2f}    Minimum {float(rec['min_score']):.2f}",
         ha="center",va="center",fontsize=9.1,color=C["sub"])
 ax.text(0,-.38,"RBH-linked rank relative to focal gene",ha="center",va="center",fontsize=8.7,color=C["sub"])
plt.subplots_adjust(left=.08,right=.985,bottom=.078,top=.885,hspace=.45,wspace=.13)
fig.suptitle("Supplementary Fig. 1  |  Local gene-order conservation across apple, pear and loquat",
             x=.076,ha="left",y=.983,fontsize=16,fontweight="bold",color=C["ink"])
fig.text(.083,.940,"Tracks show +/-5 RBH ranks; scores use +/-10 target ranks. * Local track reversed for alignment.",
         ha="left",fontsize=10,color=C["sub"])
fig.legend(handles=[Patch(facecolor=C["focal"],label="Focal locus"),
        Patch(facecolor=C["mint"],label="Apple-neighbor anchor"),
        Patch(facecolor=C["gray"],label="Other target-neighbor gene")],
        loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.53,.016),fontsize=9.5)
savefig(fig,"FigS1_V10_LocalCollinearity");plt.close(fig)
print("Wrote Fig1 and FigS1 in",OUT)
