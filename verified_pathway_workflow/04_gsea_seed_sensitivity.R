#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(clusterProfiler))
ROOT<-path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
OUT<-path.expand(Sys.getenv("LOQUAT_AUDIT_DIR",unset=file.path(ROOT,"16_PATHWAY_PERMUTATION_AUDIT_20261022")))
ttg<-read.csv(file.path(OUT,"KEGG_MDM_FULL_TERM2GENE_20261022.csv"))
tab<-read.csv(file.path(ROOT,"06_meta_analysis/apple_main36/DE_stage1_M9_vs_MM106.csv"),check.names=FALSE)
tab<-tab[is.finite(tab$stat)&!is.na(tab$GENEID),]
tab$GENEID<-as.character(tab$GENEID)
tab<-tab[!duplicated(tab$GENEID),]
gl<-sort(setNames(tab$stat,tab$GENEID),decreasing=TRUE)
# Stage2 seed already tested and robust for mdm00650; evaluate stage1 boundary at four alternate seeds.
seeds<-c(1,42,20261022,314159,271828)
l<-list()
for(seed in seeds){
 set.seed(seed)
 x<-clusterProfiler::GSEA(geneList=gl,TERM2GENE=ttg[,c("pathway","gene")],
 minGSSize=10,maxGSSize=500,pvalueCutoff=1,verbose=FALSE,seed=TRUE)
 z<-as.data.frame(x)
 hit<-z[z$ID=="mdm00650",,drop=FALSE]
 l[[length(l)+1]]<-data.frame(seed=seed,stage1_mdm00650_FDR=hit$p.adjust,
 stage1_mdm00650_NES=hit$NES,stage1_significant_pathway_n=sum(z$p.adjust<.1))
}
ans<-do.call(rbind,l)
write.csv(ans,file.path(OUT,"GSEA_SEED_SENSITIVITY.csv"),row.names=FALSE)
print(ans)
