#!/usr/bin/env Rscript
ROOT<-path.expand(Sys.getenv("LOQUAT_AUDIT_DIR",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC/16_PATHWAY_PERMUTATION_AUDIT_20261022"))
z<-readRDS(file.path(ROOT,"PERMUTATION_ES_MATRICES.RDS"))
ids147<-rownames(z$S1)
sel<-list(historical37=z$selectedOLD,frozen38=z$selectedNEW,all147=ids147)
results<-list()
for(n in names(sel)){
 ids<-sel[[n]]; f<-function(i,j)sum(sign(z$S1[ids,i])*sign(z$S2[ids,j])<0)
 obs<-f(1,1)
 cross<-unlist(lapply(1:20,function(i)vapply(1:20,function(j)f(i,j),numeric(1))))
 paired<-vapply(1:20,function(i)f(i,i),numeric(1))
 # preserving a shared biological unit across both stage measurements is appropriate ONLY if designs are paired.
 results[[n]]<-data.frame(subset=n,n_pathways=length(ids),observed=obs,cross_stage_unrestricted_perms=length(cross),
  cross_stage_unrestricted_extreme=sum(cross>=obs),cross_stage_unrestricted_P=mean(cross>=obs),
  coupled_label_perms=length(paired),coupled_label_extreme=sum(paired>=obs),
  coupled_label_P=mean(paired>=obs),mean_cross=mean(cross),mean_coupled=mean(paired))
}
r<-do.call(rbind,results)
write.csv(r,file.path(ROOT,"INDEPENDENT_VS_COUPLED_PERMUTATION.csv"),row.names=FALSE)
print(r,row.names=FALSE)
