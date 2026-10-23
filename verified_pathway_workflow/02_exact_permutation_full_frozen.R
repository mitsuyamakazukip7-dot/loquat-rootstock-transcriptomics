#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(limma);library(fgsea)})
ROOT<-path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
OUT<-path.expand(Sys.getenv("LOQUAT_AUDIT_DIR",unset=file.path(ROOT,"16_PATHWAY_PERMUTATION_AUDIT_20261022")))
O<-file.path(ROOT,"06_meta_analysis/apple_main36")
v<-read.csv(file.path(O,"vst_biological_matrix.csv"),check.names=FALSE)
genes<-as.character(v[[1]]);expr<-as.matrix(v[,-1]);rownames(expr)<-genes
ttg<-read.csv(file.path(OUT,"KEGG_MDM_FULL_TERM2GENE_20261022.csv"),stringsAsFactors=FALSE)
sets<-split(ttg$gene,ttg$pathway)
sets<-lapply(sets,intersect,y=genes)
sets<-sets[lengths(sets)>=10 & lengths(sets)<=500]
old<-read.csv(file.path(O,"cross_stage_KEGG_concordance.csv"),check.names=FALSE)
new<-read.csv(file.path(OUT,"GSEA_FULL_FROZEN_STAGE_COMPARISON.csv"),check.names=FALSE)
selectedOLD<-intersect(old$ID[old$p.adjust_s1<.1 & old$p.adjust_s2<.1],names(sets))
selectedNEW<-intersect(new$ID[new$p.adjust_s1<.1 & new$p.adjust_s2<.1],names(sets))
allids<-intersect(intersect(as.character(old$ID),as.character(new$ID)),names(sets))
stopifnot(length(sets)==147,length(selectedOLD)==37,length(selectedNEW)==38,length(allids)==147)
s1<-c("M9_S1_B1","M9_S1_B2","M9_S1_B3","MM106_S1_B1","MM106_S1_B2","MM106_S1_B3")
s2<-c("M9_S2_B1","M9_S2_B2","M9_S2_B3","MM106_S2_B1","MM106_S2_B2","MM106_S2_B3")
comb<-combn(6,3,simplify=FALSE)
stopifnot(identical(comb[[1]],1:3))
makeES<-function(mat) {
 do.call(cbind,lapply(comb,function(grpidx){
  group<-factor(ifelse(seq_len(6)%in%grpidx,"M9","MM106"),levels=c("MM106","M9"))
  fit<-eBayes(lmFit(mat,model.matrix(~group)),trend=TRUE)
  ranks<-sort(setNames(fit$t[,2],rownames(mat)),decreasing=TRUE)
  vapply(sets,function(g){
   idx<-match(g,names(ranks));idx<-idx[!is.na(idx)]
   fgsea::calcGseaStat(ranks,idx)
  },numeric(1))
 }))
}
cat("Scoring stage 1...\n")
S1<-makeES(expr[,s1,drop=FALSE])
cat("Scoring stage 2...\n")
S2<-makeES(expr[,s2,drop=FALSE])
rownames(S1)<-rownames(S2)<-names(sets)
colnames(S1)<-colnames(S2)<-sprintf("assignment%02d",seq_along(comb))
saveRDS(list(S1=S1,S2=S2,sets=sets,comb=comb,selectedOLD=selectedOLD,selectedNEW=selectedNEW),file.path(OUT,"PERMUTATION_ES_MATRICES.RDS"))
rev<-function(i,j,ids)sum(sign(S1[ids,i])*sign(S2[ids,j])<0)
rho<-function(i,j,ids)cor(S1[ids,i],S2[ids,j],method="spearman",use="complete.obs")
pairs<-expand.grid(i=seq_along(comb),j=seq_along(comb))
result<-list()
for(label in c("historical_selected37","frozen_full_selected38","all147")){
 ids<-switch(label,historical_selected37=selectedOLD,frozen_full_selected38=selectedNEW,all147=allids)
 cnt<-vapply(seq_len(nrow(pairs)),function(k)rev(pairs$i[k],pairs$j[k],ids),integer(1))
 obs<-rev(1,1,ids)
 q<-rho(1,1,ids)
 nullrho<-vapply(seq_len(nrow(pairs)),function(k)rho(pairs$i[k],pairs$j[k],ids),numeric(1))
 result[[label]]<-data.frame(subset=label,n_pathways=length(ids),observed_reversals=obs,
  n_permutation=length(cnt),tail_ge=sum(cnt>=obs),p_ge=mean(cnt>=obs),
  observed_rho=q,p_rho_le=mean(nullrho<=q),mean_null_rev=mean(cnt),sd_null_rev=sd(cnt), stringsAsFactors=FALSE)
 if(label=="historical_selected37"){
  comparison<-data.frame(KEGG_pathway=ids,observed_ES_s1=S1[ids,1],observed_ES_s2=S2[ids,1],
  observed_reversal=sign(S1[ids,1])!=sign(S2[ids,1]),hist_NES_s1=old$NES_s1[match(ids,old$ID)],
  hist_NES_s2=old$NES_s2[match(ids,old$ID)],historical_reversal=sign(old$NES_s1[match(ids,old$ID)])!=sign(old$NES_s2[match(ids,old$ID)]))
  write.csv(comparison,file.path(OUT,"HISTORICAL37_GSEA_vs_LIMMA_ES.csv"),row.names=FALSE)
 }
 if(label=="all147"){
  oldrho<-read.csv(file.path(O,"../reviewer_v8/pathway_label_permutation_summary.csv"))
  cat("Old result historical all 147 P =",oldrho$exact_permutation_p[2],"\n")
 }
 write.csv(data.frame(group=label,pair_i=pairs$i,pair_j=pairs$j,reverse_count=cnt,spearman_rho=nullrho),
  file.path(OUT,paste0("NULL400_",label,".csv")),row.names=FALSE)
}
o<-do.call(rbind,result)
write.csv(o,file.path(OUT,"PERMUTATION_REPRODUCTION_COMPARISON.csv"),row.names=FALSE)
print(o,row.names=FALSE)
cat("HISTORICAL 37 (GSEA) REVERSED",sum(sign(old$NES_s1[match(selectedOLD,old$ID)])!=sign(old$NES_s2[match(selectedOLD,old$ID)])),"\n")
