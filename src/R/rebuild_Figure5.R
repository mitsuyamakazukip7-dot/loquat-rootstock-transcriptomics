#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid);library(scales);library(svglite)})
root <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
path <- file.path(root,"12_submission_v10_assets")
out_dir <- path.expand(Sys.getenv("LOQUAT_FIGURE_OUT",unset="./output"))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
ink<-"#263238";muted<-"#66717A";lt<-"#E9ECEB";teal<-"#8BD1C8";pink<-"#F2BCD4"
gold<-"#D8BC72";tealDark<-"#438D85";pinkDark<-"#B86E8F";border<-"#CFD8D7"
base_theme<-theme_minimal(base_family="Helvetica",base_size=10)+
 theme(plot.background=element_rect(fill="white",color=NA),panel.background=element_rect(fill="white",color=NA),
 panel.grid=element_blank(),axis.text=element_text(color=ink),axis.title=element_blank(),
 plot.title=element_text(face="bold",color=ink,size=13,margin=margin(b=8)),
 plot.subtitle=element_text(color=muted,size=9,margin=margin(b=8)),
 plot.margin=margin(8,8,8,8))

ora<-read.csv(file.path(path,"Table_S5_KEGG_fixed151_all_sets.csv"),stringsAsFactors=FALSE)
fl<-ora[ora$KEGG_ID=="mdm00941",]
fl$set<-factor(fl$Candidate_Set,levels=c("S4_5","S3_24","S2_190"))
fl$desc<-factor(fl$Candidate_Set,levels=c("S4_5","S3_24","S2_190"),
                labels=c("Five genes (both tissues)","24 externally supported","190 collinearity-filtered"))
fl$score <- -log10(fl$BH_FDR_151_tested)
fl$sig<-ifelse(fl$BH_FDR_151_tested<.05,"FDR < 0.05","FDR >= 0.05")
fl$txt<-sprintf("FDR %.3g  |  %s hits / %s mapped",fl$BH_FDR_151_tested,
                fl$Genes_in_Pathway,fl$Candidates_KEGG_mapped)
threshold <- -log10(.05)
f5a<-ggplot(fl,aes(x=score,y=desc))+
 geom_vline(xintercept=threshold,linetype="dashed",color="#A8B6B4",linewidth=.5)+
 geom_segment(aes(x=0,xend=score,yend=desc),linewidth=1,color=border)+
 geom_point(aes(fill=sig),size=5.3,shape=21,color=ink,stroke=.35)+
 geom_text(aes(label=txt),hjust=0,size=3.15,nudge_x=.15,color=ink)+
 scale_fill_manual(values=c("FDR < 0.05"=teal,"FDR >= 0.05"="#C9CFCE"))+
 scale_x_continuous(limits=c(-.15,6.9),breaks=seq(0,6,1),expand=c(0,0))+
 labs(title="a  Flavonoid biosynthesis across candidate sets",
      subtitle="One-sided hypergeometric tests; BH correction across the same 151 KEGG pathways",
      x=expression(-log[10]("adjusted P value")),y=NULL)+
 base_theme+
 theme(legend.position="none",axis.text.y=element_text(size=9.7),
       plot.margin=margin(16,20,8,18))
# Panel b five-locus leave-one-out vs twenty-four-locus controls.
z<-read.csv(file.path(path,"Table_S6_PKS_paralog_sensitivity.csv"),stringsAsFactors=FALSE)
z<-z[z$Candidate_Set %in% c("S3_24","S4_5") &
       z$PKS5_Exclusion %in% c("none","PKS5","PKS5like"),]
z$panel<-factor(z$Candidate_Set,levels=c("S3_24","S4_5"),
                labels=c("24 externally supported","Five genes (both tissues)"))
z$exclude<-factor(z$PKS5_Exclusion,
                  levels=c("PKS5like","PKS5","none"),
                  labels=c("Remove PKS5-like","Remove PKS5","Both retained"))
z$score<--log10(z$BH_FDR_151_tested)
z$sig<-z$BH_FDR_151_tested<.05
z$label<-sprintf("FDR %.3g   (%s hits)",z$BH_FDR_151_tested,z$Flavonoid_Gene_Hits)
f5b<-ggplot(z,aes(x=score,y=exclude))+
 geom_vline(xintercept=threshold,linetype="dashed",color="#A8B6B4",linewidth=.5)+
 geom_segment(aes(x=0,xend=score,yend=exclude),color=border,linewidth=.9)+
 geom_point(aes(fill=sig),size=4.5,shape=21,color=ink,stroke=.3)+
 geom_text(aes(label=label),hjust=0,nudge_x=.12,size=3.05,color=ink)+
 facet_wrap(~panel,ncol=1,scales="free_y")+
 scale_fill_manual(values=c("TRUE"=teal,"FALSE"="#C9CFCE"))+
 scale_x_continuous(limits=c(-.12,4.4),breaks=0:4,expand=c(0,0))+
 labs(title="b  Dependence on adjacent PKS5-related loci",
      subtitle="Removing either gene reduces the flavonoid hit count from two to one",
      x=expression(-log[10]("adjusted P value")),y=NULL)+
 base_theme+
 theme(legend.position="none",strip.text=element_text(face="bold",size=9),
       axis.text.y=element_text(size=8.6),
       plot.margin=margin(16,15,8,15))
f5<-f5a/f5b+plot_layout(heights=c(2.25,3.4))
ggsave(file.path(out_dir,"Figure5_Rebuilt.png"),f5,
       width=11.2,height=9.6,dpi=300,bg="white")
ggsave(file.path(out_dir,"Figure5_Rebuilt.pdf"),f5,
       width=11.2,height=9.6,device="pdf",bg="white")
ggsave(file.path(out_dir,"Figure5_Rebuilt.svg"),f5,
       width=11.2,height=9.6,device=svglite,bg="white")
message("Figure 5 rendered to ",out_dir)
