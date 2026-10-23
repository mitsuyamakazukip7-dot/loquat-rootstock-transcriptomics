#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
root <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
source <- path.expand(Sys.getenv("LOQUAT_FIG3_SOURCE",unset=file.path(root,"12_submission_v10_assets","Figure3_Original_9Pathway_Source.csv")))
output <- path.expand(Sys.getenv("LOQUAT_FIGURE_OUT",unset="./output"))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(source,stringsAsFactors=FALSE)
stopifnot(nrow(d)==9)
call <- function(nes,q) {if(is.na(nes) || is.na(q)) return("not measured")
 if(q>=0.1) return("ns")
 if(nes>0) return("Up")
 if(nes<0) return("Down")
 return("ns")}
r <- list();k <- 1
for(i in seq_len(nrow(d))){
 a<-d[i,]
 for(j in 1:3){
  t<-if(j==1)ifelse(a$pear_Q<.05,ifelse(a$pear_Z>0,"Up","Down"),"ns") else if(j==2)call(a$apple_s1_NES,a$apple_s1_FDR) else call(a$apple_s2_NES,a$apple_s2_FDR)
  r[[k]]<-data.frame(pathway=a$pathway,experiment=c("Pear","Apple stage 1","Apple stage 2")[j],status=t)
  k<-k+1
 }
}
z<-do.call(rbind,r)
z$pathway<-factor(z$pathway,levels=rev(d$pathway))
z$experiment<-factor(z$experiment,levels=c("Pear","Apple stage 1","Apple stage 2"))
p<-ggplot(z,aes(experiment,pathway,fill=status))+
 geom_tile(color="white",linewidth=0.65)+
 geom_text(aes(label=ifelse(status%in%c("not measured","ns"),ifelse(status=="ns","ns","NA"),status)),size=3.25,color="#263238")+
 scale_fill_manual(values=c("Up"="#FFF4B3","Down"="#B3DDD2","ns"="#ECEEEE","not measured"="#FFFFFF"),drop=FALSE)+
 scale_x_discrete(position="top")+
 labs(title="Pear-associated pathways and apple rootstock response contexts",
 subtitle="Pear Q < 0.05; apple GSEA FDR < 0.1. Directions and significance only.",x=NULL,y=NULL)+
 theme_minimal(base_family="Helvetica",base_size=11)+
 theme(panel.grid=element_blank(),axis.text.x=element_text(face="bold"),legend.position="bottom",
 plot.title=element_text(face="bold"),plot.margin=margin(20,20,20,20))
ggsave(file.path(output,"Figure3_From_Historical9Pathway_Source.pdf"),p,width=11,height=7,device="pdf",bg="white")
ggsave(file.path(output,"Figure3_From_Historical9Pathway_Source.svg"),p,width=11,height=7,device="svg",bg="white")
ggsave(file.path(output,"Figure3_From_Historical9Pathway_Source.png"),p,width=11,height=7,dpi=600,bg="white")
message("FIG3_VALIDATED_SOURCE_ROWS:",nrow(d))
