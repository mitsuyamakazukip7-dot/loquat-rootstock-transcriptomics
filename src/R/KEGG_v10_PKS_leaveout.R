#!/usr/bin/env Rscript
# Frozen KEGG memberships, all 151 tests; full vs leave-one-out PKS paralog controls.
P<-file.path(path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC")),"11_methods_audit_20261020")
S<-read.csv(file.path(P,"KEGG_frozen_TERM2GENE.csv"),stringsAsFactors=FALSE)
G<-split(S$gene,S$pathway);U<-scan(file.path(P,"KEGG_frozen_mapped_background.txt"),what="",quiet=TRUE)
pk<-c("103421794","103443512")
out<-list()
for(k in c("pool","external","both")){
 f<-file.path(P,"v10_universal_synteny",paste0("s0.8_b0_",k,".csv"))
 d<-read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
 for(drop in c("none","PKS5","PKS5like","both")){
  x<-as.character(d$apple_geneid)
  if(drop=="PKS5")x<-setdiff(x,pk[1])
  if(drop=="PKS5like")x<-setdiff(x,pk[2])
  if(drop=="both")x<-setdiff(x,pk)
  ids<-intersect(x,U);n<-length(ids);N<-length(U)
  q<-do.call(rbind,lapply(names(G),function(path){
    K<-length(G[[path]]);nHit<-sum(ids%in%G[[path]])
    data.frame(group=k,leaveout=drop,pathway=path,raw_set_size=length(x),
        n_mapped=n,nHit=nHit,p=if(nHit==0)1 else phyper(nHit-1,K,N-K,n,lower.tail=FALSE))
  }))
  q$FDR_BH_151<-p.adjust(q$p,method="BH")
  out[[paste(k,drop)]]<-q
 }
}
r<-do.call(rbind,out)
write.csv(r,file.path(P,"KEGG_v10_PKS_leaveout_ALL_paths.csv"),row.names=FALSE)
fl<-r[r$pathway=="mdm00941",]
write.csv(fl,file.path(P,"KEGG_v10_PKS_leaveout_summary.csv"),row.names=FALSE)
print(fl[,c("group","leaveout","raw_set_size","n_mapped","nHit","p","FDR_BH_151")],row.names=FALSE)
