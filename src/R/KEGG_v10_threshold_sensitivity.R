#!/usr/bin/env Rscript
# Re-evaluate ALL v10 threshold candidates with fixed KEGG TERM2GENE.
# No annotation membership is used in candidate ranking or filtering.
ROOT <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
base <- file.path(ROOT,"11_methods_audit_20261020")
input <- file.path(base,"v10_universal_synteny")
ttg <- read.csv(file.path(base,"KEGG_frozen_TERM2GENE.csv"),stringsAsFactors=FALSE)
stopifnot(all(c("pathway","gene") %in% names(ttg)))
universe <- scan(file.path(base,"KEGG_frozen_mapped_background.txt"),what="",quiet=TRUE)
sets <- lapply(split(as.character(ttg$gene),as.character(ttg$pathway)),unique)
all <- read.csv(file.path(input,"eligible_269.csv"),check.names=FALSE)
legacy <- read.csv(file.path(ROOT,"09_tables/Table_candidate_blind_ranking_all20795.csv"),check.names=FALSE)
evidence_score <- setNames(legacy$blind_base_score,legacy$apple_geneid)
all$blind_base_score_from_saved <- as.numeric(evidence_score[as.character(all$apple_geneid)])
stopifnot(nrow(all)==269,length(sets)==151,length(universe)==4775)
def <- list()
for(th in c(.6,.7,.8,.9,1.0)){
  z <- all[all$n_reference_apple_anchors>=8 & all$synteny_score_v10_min>=th,]
  for(minpoints in c(0,11)){
   v <- z[!is.na(z$blind_base_score_from_saved) & z$blind_base_score_from_saved>=minpoints,]
   for(group in c("pool","external","both")){
     g <- if(group=="pool")v else if(group=="external")v[tolower(as.character(v$external_any))=="true",] else v[tolower(as.character(v$external_both))=="true",]
     key <- sprintf("s%.1f_b%d_%s",th,minpoints,group)
     write.csv(g,file.path(input,paste0(key,".csv")),row.names=FALSE)
     def[[key]] <- as.character(g$apple_geneid)
   }
  }
}
results<-list()
for(key in names(def)){
 ids <- intersect(unique(def[[key]]),universe)
 N<-length(universe);n<-length(ids)
 sub<-do.call(rbind,lapply(names(sets),function(path){
   K<-length(sets[[path]]);k<-sum(ids%in%sets[[path]])
   p<-if(n==0 || k==0)1 else phyper(k-1,K,N-K,n,lower.tail=FALSE)
   data.frame(id=key,pathway=path,original_n=length(def[[key]]),mapped_n=n,hits=k,raw_p=p,
     fold=ifelse(k>0 && n>0,(k/n)/(K/N),0),stringsAsFactors=FALSE)
 }))
 sub$FDR_BH_151<-p.adjust(sub$raw_p,method="BH")
 results[[key]]<-sub
}
dat <- do.call(rbind,results)
write.csv(dat,file.path(base,"KEGG_v10_synteny_full_FDR.csv"),row.names=FALSE)
fl<-dat[dat$pathway=="mdm00941",]
fl<-fl[order(fl$id),]
write.csv(fl,file.path(base,"KEGG_v10_flavonoid_summary.csv"),row.names=FALSE)
print(fl[,c("id","original_n","mapped_n","hits","fold","FDR_BH_151")],row.names=FALSE)
