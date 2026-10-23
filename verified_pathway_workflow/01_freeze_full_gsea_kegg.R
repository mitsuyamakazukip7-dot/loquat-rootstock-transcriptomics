#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(KEGGREST);library(clusterProfiler)})
ROOT<-path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
OUT<-path.expand(Sys.getenv("LOQUAT_AUDIT_DIR",unset=file.path(ROOT,"16_PATHWAY_PERMUTATION_AUDIT_20261022")))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
TERM<-file.path(OUT,"KEGG_MDM_FULL_TERM2GENE_20261022.csv")
if(!file.exists(TERM)){
 links<-KEGGREST::keggLink("pathway","mdm")
 table<-unique(data.frame(pathway=sub("^path:","",unname(links)),gene=sub("^mdm:","",names(links)),stringsAsFactors=FALSE))
 table<-table[order(table$pathway,table$gene),]
 write.csv(table,TERM,row.names=FALSE)
}
ttg<-read.csv(TERM,stringsAsFactors=FALSE)
stopifnot(length(unique(ttg$pathway))>=150)
v<-read.csv(file.path(ROOT,"06_meta_analysis/apple_main36/vst_biological_matrix.csv"),check.names=FALSE)
genes<-as.character(v[[1]])
term<-split(ttg$gene,ttg$pathway)
lens<-sapply(term,function(g)length(intersect(g,genes)))
status<-data.frame(pathway=names(lens),raw_Kegg_n=lengths(term),measured_n=as.numeric(lens),eligible10to500=lens>=10&lens<=500)
write.csv(status,file.path(OUT,"KEGG_SET_SIZE_BY_UNIVERSE.csv"),row.names=FALSE)
print(c(links=nrow(ttg),pathways=length(term),full_measured_10to500=sum(lens>=10 & lens<=500)))
# pin original tested lists, to distinguish re-evaluation of selected pathways from re-selection
old1<-read.csv(file.path(ROOT,"06_meta_analysis/apple_main36/stage1_gseKEGG.csv"),check.names=FALSE)
old2<-read.csv(file.path(ROOT,"06_meta_analysis/apple_main36/stage2_gseKEGG.csv"),check.names=FALSE)
term2name<-data.frame(pathway=old1$ID,name=old1$Description)
write.csv(term2name,file.path(OUT,"KEGG_TERM2NAME_HISTORICAL_LABELS.csv"),row.names=FALSE)
run_stage<-function(stage, seed) {
 raw<-read.csv(file.path(ROOT,"06_meta_analysis/apple_main36",sprintf("DE_stage%d_M9_vs_MM106.csv",stage)),check.names=FALSE)
 y<-raw[is.finite(raw$stat)&!is.na(raw$GENEID),]
 y$GENEID<-as.character(y$GENEID)
 y<-y[!duplicated(y$GENEID),]
 gl<-sort(setNames(y$stat,y$GENEID),decreasing=TRUE)
 set.seed(seed)
 fit<-clusterProfiler::GSEA(geneList=gl,TERM2GENE=ttg[,c("pathway","gene")],
   TERM2NAME=term2name[,c("pathway","name")],minGSSize=10,maxGSSize=500,pvalueCutoff=1,
   verbose=FALSE,seed=TRUE,by="fgsea")
 rr<-as.data.frame(fit)
 write.csv(rr,file.path(OUT,sprintf("GSEA_FULL_FROZEN_STAGE%d.csv",stage)),row.names=FALSE)
 rr
}
a<-run_stage(1,20261022)
b<-run_stage(2,20261023)
x<-merge(a[,c("ID","setSize","NES","p.adjust")],b[,c("ID","setSize","NES","p.adjust")],by="ID",suffixes=c("_s1","_s2"))
x$sharedFDR10<-x$p.adjust_s1<.1 & x$p.adjust_s2<.1
x$reversed<-sign(x$NES_s1)!=sign(x$NES_s2)
write.csv(x,file.path(OUT,"GSEA_FULL_FROZEN_STAGE_COMPARISON.csv"),row.names=FALSE)
old<-merge(old1[,c("ID","setSize","NES","p.adjust")],old2[,c("ID","setSize","NES","p.adjust")],by="ID",suffixes=c("_s1","_s2"))
old$sharedFDR10<-old$p.adjust_s1<.1 & old$p.adjust_s2<.1
old$reversed<-sign(old$NES_s1)!=sign(old$NES_s2)
o<-merge(old,x,by="ID",all=TRUE,suffixes=c("_old","_new"))
write.csv(o,file.path(OUT,"GSEA_OLD_VS_FROZEN_FULL_COMPARISON.csv"),row.names=FALSE)
cat("OLD total",nrow(old),"shared",sum(old$sharedFDR10),"reversed",sum(old$sharedFDR10&old$reversed),"\n")
cat("NEW total",nrow(x),"shared",sum(x$sharedFDR10),"reversed",sum(x$sharedFDR10&x$reversed),"\n")
cat("SET SIZE OLD VS NEW FOR FLAVONOID",old$setSize_s1[old$ID=="mdm00941"],x$setSize_s1[x$ID=="mdm00941"],"\n")
cat("HIGH CONFIDENCE old/new old-selected-overlap",sum(o$sharedFDR10_old&o$sharedFDR10_new,na.rm=TRUE),"\n")
