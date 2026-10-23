#!/usr/bin/env Rscript
suppressPackageStartupMessages({
library(ggplot2);library(patchwork);library(grid);library(scales);library(svglite)
})
root <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
path <- file.path(root,"12_submission_v10_assets")
out_dir <- path.expand(Sys.getenv("LOQUAT_FIGURE_OUT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC/16_FIGURE_REVISED_20261022")) 
d <- read.csv(file.path(path,"Figure4_24_heatmap_source.csv"),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(d)==24,
          all(tolower(as.character(d$External_any_FDR_lt_0p05))=="true"))
# Colors match the manuscript's restrained pastel palette.
ink<-"#263238";muted<-"#66717A";lt<-"#E9ECEB";teal<-"#8BD1C8";pink<-"#F2BCD4"
gold<-"#D8BC72";tealDark<-"#438D85";pinkDark<-"#B86E8F";border<-"#CFD8D7"
base_theme<-theme_minimal(base_family="Helvetica",base_size=10)+
 theme(plot.background=element_rect(fill="white",color=NA),
       panel.background=element_rect(fill="white",color=NA),
       panel.grid=element_blank(),
       axis.text=element_text(color=ink),axis.title=element_blank(),
       plot.title=element_text(face="bold",color=ink,size=13,margin=margin(b=8)),
       plot.subtitle=element_text(color=muted,size=9,margin=margin(b=8)),
       plot.margin=margin(8,8,8,8))
# Fig4a: flow, annotation-free filters; annotations are inspected AFTER selection.
stage_names<-c("Apple-centered RBH","Apple-pear direction +\nloquat stem expression",
"Local collinearity\nin both species","External apple\nat least one tissue",
"External apple\nboth tissues")
nums<-c("20,795","269","190","24","5")
subtxt<-c("One-to-one correspondences","Both expression criteria","v10 score >= 0.80","FDR < 0.05","Same sign + FDR < 0.05")
flow<-data.frame(i=1:5,x=seq(0.6,9.4,length.out=5),label=stage_names,num=nums,sub=subtxt)
flow$xmin<-flow$x-.86;flow$xmax<-flow$x+.86
f4a<-ggplot(flow)+
 geom_segment(data=data.frame(x=flow$x[-5]+.9,xend=flow$x[-1]-.9),
              aes(x=x,xend=xend,y=0,yend=0),
              linewidth=.6,color="#91A6A3",
              arrow=arrow(length=unit(2.2,"mm"),type="closed"))+
 geom_rect(aes(xmin=xmin,xmax=xmax,ymin=-.60,ymax=.60),
           fill=c("#F0F3F3","#E1EFED","#B3DDD2","#8BD1C8","#D8BC72"),
           color=border,linewidth=.5)+
 geom_text(aes(x=x,y=.25,label=num),size=5,fontface="bold",color=ink)+
 geom_text(aes(x=x,y=-.05,label=label),size=2.6,color=ink,lineheight=.98)+
 geom_text(aes(x=x,y=-.43,label=sub),size=2.15,color=muted)+
 coord_cartesian(xlim=c(-.4,10.4),ylim=c(-.7,.85),expand=FALSE,clip="off")+
 annotate("text",x=-.3,y=.85,label="a",fontface="bold",size=5.2,hjust=0,vjust=1,color=ink)+
 theme_void(base_family="Helvetica")+
 theme(plot.margin=margin(18,22,0,16),plot.background=element_rect(fill="white",color=NA))
# Fig4b: expression effects are colored ONLY if meeting the appropriate support criterion.
# Pear = direction + large-effect selection, NOT pear biological-replicate inference.
d$gene_short <- d$Apple_Annotation_postselection
d$gene_short <- gsub(", transcript variant.*","",d$gene_short)
d$gene_short <- gsub("probable ","",d$gene_short,fixed=TRUE)
d$gene_short <- gsub("protein ","",d$gene_short,fixed=TRUE)
d$gene_short <- gsub("polyketide synthase 5-like","PKS5-like",d$gene_short,fixed=TRUE)
d$gene_short <- gsub("polyketide synthase 5","PKS5",d$gene_short,fixed=TRUE)
d$gene_short <- gsub("galactinol--sucrose galactosyltransferase 1","galactinol-sucrose GT1",d$gene_short,fixed=TRUE)
d$gene_short <- gsub("chloroplastic","chloropl.",d$gene_short,fixed=TRUE)
d$gene_short <- substr(d$gene_short,1,38)
# Row labels disambiguate similarly named genes.
d$gene_label <- paste0(d$gene_short,"  [",d$Apple_GeneID,"]")
d$gene_label <- make.unique(d$gene_label)
# keep both-tissue subgroup at top
d$gene_label <- factor(d$gene_label,levels=rev(d$gene_label))
make_stat<-function(fdr,fc,pvalueCut=.05,no_q_filter=FALSE){
 if(!is.finite(fc))return("Not significant")
 if(no_q_filter)return(ifelse(fc>0,"Up","Down"))
 if(is.finite(fdr) && fdr<pvalueCut) return(ifelse(fc>0,"Up","Down"))
 return("Not significant")
}
labs<-c("Apple S1","Apple S2","Pear effect","Ext. new-shoot\nleaves","Ext. fruiting-spur\nleaves")
rows<-list(); i<-1
for(k in seq_len(nrow(d))){
 r<-d[k,]
 vals<-list(c(r$Apple_Stage1_FDR,r$Apple_Stage1_log2FC),
            c(r$Apple_Stage2_FDR,r$Apple_Stage2_log2FC),
            c(NA,r$Pear_log2FC),
            c(r$External_shoot_FDR,r$External_shoot_log2FC),
            c(r$External_spur_FDR,r$External_spur_log2FC))
 for(j in seq_along(vals)){
   v<-vals[[j]]
   rows[[i]]<-data.frame(gene_label=r$gene_label,assay=labs[j],
                        class=make_stat(v[1],v[2],no_q_filter=(j==3)),
                        n=k);i<-i+1
 }
}
tiles<-do.call(rbind,rows)
tiles$assay<-factor(tiles$assay,levels=labs)
tiles$class<-factor(tiles$class,levels=c("Down","Not significant","Up"))
side<-data.frame(gene_label=d$gene_label,
      x=rep(6.06,nrow(d)),syn=sprintf("%.2f",d$V10_min_collinearity_score),
      tpm=sprintf("%.1f",d$Loquat_Stem_median_TPM),
      subgroup=ifelse(d$External_both_same_direction_FDR_lt_0p05,"★",""))
side$y<-as.numeric(d$gene_label)
# manually construct one combined ggplot with five colored tests and two text columns.
tiles$xn<-as.numeric(tiles$assay)
f4b<-ggplot()+
 geom_tile(data=tiles,aes(x=xn,y=gene_label,fill=class),height=.83,width=.80,
           color="#FFFFFF",linewidth=.40)+
 geom_text(data=side,aes(x=6.04,y=gene_label,label=syn),size=2.65,color=ink)+
 geom_text(data=side,aes(x=7.12,y=gene_label,label=tpm),size=2.65,color=ink)+
 geom_point(data=side[side$subgroup=="★",],aes(x=.45,y=gene_label),
           color=tealDark,size=1.6,shape=16)+
 annotate("point",x=.43,y=25.1,size=1.6,color=tealDark)+
 scale_fill_manual(values=c("Down"=pink,"Not significant"="#E8ECEC","Up"=teal),
                   drop=FALSE)+
 scale_x_continuous(limits=c(.1,7.65),breaks=c(1:5,6.04,7.12),
                    labels=c(labs,"Min. synteny","Stem TPM"),position="top",expand=c(0,0))+
 scale_y_discrete(expand=expansion(add=c(.45,1.15)))+
 guides(fill=guide_legend(override.aes=list(color=NA),title="Effect direction"))+
 labs(title="b  Evidence matrix for the 24 externally supported loci",
      subtitle="Expression colors indicate the observed effect direction; gray indicates FDR >= 0.05. Pear colors indicate direction with |log2FC| >= 1 (no independent biological replication).")+
 base_theme+
 theme(axis.text.x=element_text(size=8.7,face="bold",hjust=.5,margin=margin(b=5)),
       axis.text.y=element_text(size=7.75,color=ink),
       legend.position="bottom",legend.title=element_text(size=8),
       legend.text=element_text(size=8),
       plot.title=element_text(size=12.5,face="bold"),
       plot.subtitle=element_text(size=8.3),
       plot.margin=margin(2,13,8,12))
f4<-f4a/f4b+plot_layout(heights=c(1.2,6.8))
ggsave(file.path(out_dir,"Figure4_REVISED.png"),f4,
       width=14.0,height=10.6,dpi=600,bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"Figure4_REVISED.pdf"),f4,
       width=14.0,height=10.6,device="pdf",bg="white",limitsize=FALSE)
ggsave(file.path(out_dir,"Figure4_REVISED.svg"),f4,
       width=14.0,height=10.6,device=svglite,bg="white",limitsize=FALSE)


message("FIG4_RENDER_SUCCESS: 24 matrix rows, 269 -> 190 -> 24 -> 5")
