#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(scales); library(grid)
})

ROOT <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
OUT <- path.expand(Sys.getenv("LOQUAT_FIGURE_OUT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC/16_FIGURE_REVISED_20261022"))
PAN <- file.path(OUT, "panels"); COMB <- file.path(OUT, "combined")
dir.create(PAN, recursive=TRUE, showWarnings=FALSE)
dir.create(COMB, recursive=TRUE, showWarnings=FALSE)

# Requested palette
C_MINT_LIGHT <- "#b3ddd2"
C_MINT       <- "#8bd1c8"
C_PINK       <- "#f2bcd4"
C_YELLOW     <- "#fff4b3"
C_GOLD       <- "#d8bc72"
C_TEXT       <- "#263238"
C_GRAY       <- "#6b7280"
C_LIGHT      <- "#e5e7eb"
C_BG         <- "#fafafa"

nature_theme <- theme_classic(base_size=7, base_family="Helvetica") +
  theme(axis.line=element_line(linewidth=0.35, colour=C_TEXT),
        axis.ticks=element_line(linewidth=0.35, colour=C_TEXT),
        axis.text=element_text(colour=C_TEXT, size=6.2),
        axis.title=element_text(colour=C_TEXT, size=7),
        legend.text=element_text(size=5.8),
        legend.title=element_blank(),
        plot.title=element_text(size=8, face="bold", colour=C_TEXT, hjust=0),
        plot.subtitle=element_text(size=6.2, colour=C_GRAY),
        plot.margin=margin(5,5,5,5))

tag <- function(p, x) p + labs(tag=x) + theme(plot.tag=element_text(face="bold", size=9), plot.tag.position=c(0.01,0.99))

save_both <- function(p, stem, w_mm, h_mm, dir=PAN) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width=w_mm/25.4, height=h_mm/25.4, units="in", device="pdf", bg="white")
  ggsave(file.path(dir, paste0(stem, ".png")), p, width=w_mm/25.4, height=h_mm/25.4, units="in", dpi=600, bg="white")
}

# ---------------- Fig 1 ----------------
blank <- theme_void(base_family="Helvetica") + theme(plot.margin=margin(5,5,5,5))

fig1a <- ggplot() +
  # vigorous tree
  annotate("segment", x=.20,xend=.20,y=.20,yend=.58, linewidth=2.0, colour="#8b6a4b", lineend="round") +
  annotate("segment", x=.20,xend=.10,y=.38,yend=.52, linewidth=.65, colour="#8b6a4b") +
  annotate("segment", x=.20,xend=.30,y=.38,yend=.52, linewidth=.65, colour="#8b6a4b") +
  annotate("point", x=c(.10,.16,.20,.24,.30,.14,.26), y=c(.54,.63,.67,.63,.54,.72,.72), size=c(12,13,14,13,12,11,11), colour=C_MINT_LIGHT) +
  # compact tree
  annotate("segment", x=.80,xend=.80,y=.20,yend=.47, linewidth=1.7, colour="#8b6a4b", lineend="round") +
  annotate("segment", x=.80,xend=.73,y=.33,yend=.43, linewidth=.55, colour="#8b6a4b") +
  annotate("segment", x=.80,xend=.87,y=.33,yend=.43, linewidth=.55, colour="#8b6a4b") +
  annotate("point", x=c(.73,.78,.82,.87,.76,.84), y=c(.46,.53,.53,.46,.59,.59), size=c(9,10,10,9,8,8), colour=C_MINT) +
  # graft/root indication
  annotate("segment", x=.80,xend=.80,y=.11,yend=.21, linewidth=2.1, colour=C_GOLD, lineend="round") +
  annotate("point", x=.80,y=.20,size=3.4,shape=21,fill="white",colour=C_GOLD,stroke=.8) +
  annotate("segment", x=.50,xend=.50,y=.53,yend=.66, linewidth=2.0, colour="#8b6a4b") +
  annotate("segment", x=.50,xend=.50,y=.35,yend=.45, linewidth=2.0, colour=C_GOLD) +
  annotate("text",x=.50,y=.49,label="+",size=5.5,colour=C_GRAY) +
  annotate("segment",x=.58,xend=.69,y=.49,yend=.49,arrow=arrow(length=unit(2,"mm")),linewidth=.45,colour=C_GOLD) +
  annotate("text", x=.20,y=.09,label="Vigorous / standard",size=2.5,fontface="bold",colour=C_TEXT) +
  annotate("text", x=.80,y=.09,label="Dwarfing rootstock",size=2.5,fontface="bold",colour=C_GOLD) +
  annotate("text", x=.50,y=.70,label="same fruiting scion",size=2.35,fontface="bold",colour=C_TEXT) +
  annotate("text", x=.50,y=.31,label="dwarfing rootstock",size=2.35,fontface="bold",colour=C_GOLD) +
  annotate("text", x=.84,y=.22,label="graft union",size=2.2,hjust=0,colour=C_GOLD) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.82),clip="off") +
  labs(title="Dwarfing changes tree architecture while retaining the scion") + blank +
  theme(plot.title=element_text(size=8,face="bold",colour=C_TEXT,hjust=0))
fig1a <- tag(fig1a,"a")

benefits <- tibble(x=c(.12,.37,.62,.87), label=c("Compact\ncanopy","High-density\nplanting","Easier pruning\n& harvest","Manageable orchard\narchitecture"), fill=c(C_MINT_LIGHT,C_MINT,C_PINK,C_YELLOW), edge=c(C_MINT,C_MINT,C_PINK,C_GOLD))
fig1b <- ggplot(benefits) +
  geom_tile(aes(x=x,y=.50,fill=I(fill),colour=I(edge)),width=.20,height=.47,linewidth=.55) +
  geom_point(aes(x=x,y=.61,colour=I(edge)),shape=21,fill="white",size=6,stroke=.7) +
  geom_text(aes(x=x,y=.61,label="+",colour=I(edge)),fontface="bold",size=4.2) +
  geom_text(aes(x=x,y=.43,label=label),size=2.45,fontface="bold",lineheight=.95,colour=C_TEXT) +
  annotate("text",x=.5,y=.14,label="Architecture-focused benefits motivate discovery of conserved vigor-associated genes.",size=2.25,colour=C_GRAY) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.82),clip="off") +
  labs(title="Why dwarfing matters in orchard systems") + blank +
  theme(plot.title=element_text(size=8,face="bold",colour=C_TEXT,hjust=0))
fig1b <- tag(fig1b,"b")

nodes <- tibble(x=c(.16,.50,.84,.50,.27,.72,.80), y=c(.72,.72,.72,.45,.20,.20,.43),
                lab=c("Apple primary\nM9 vs MM106\n2 phenological contexts","Pear support\ndwarf vs Cuiguan\ndirectional agreement","Loquat anchor\nstem expression\nprojection target","20,795 strict 1:1:1 RBH triples\nsequence mapping + direction + expression + synteny","38-gene annotation-blind\ncandidate pool","17 external-apple-supported\nhigh-confidence genes","External apple\nindependent support"),
                fill=c("#eef9f7",C_YELLOW,"#fdf1f7",C_BG,"#f1faf8","#fff8df",C_YELLOW),
                edge=c(C_MINT,C_GOLD,C_PINK,C_LIGHT,C_MINT,C_GOLD,C_GOLD))
fig1c <- ggplot() +
  geom_label(data=nodes[1:3,],aes(x=x,y=y,label=lab,fill=I(fill)),colour=C_TEXT,size=2.25,linewidth=.45,label.r=unit(2,"mm"),fontface="bold",lineheight=.95) +
  geom_label(data=nodes[4:7,],aes(x=x,y=y,label=lab,fill=I(fill)),colour=C_TEXT,size=2.15,linewidth=.45,label.r=unit(2,"mm"),lineheight=.95) +
  annotate("segment",x=c(.16,.50,.84),xend=.50,y=c(.62,.62,.62),yend=.52,arrow=arrow(length=unit(1.7,"mm")),linewidth=.4,colour=c(C_MINT,C_GOLD,C_PINK)) +
  annotate("segment",x=.50,xend=.27,y=.37,yend=.27,arrow=arrow(length=unit(1.7,"mm")),linewidth=.4,colour=C_MINT) +
  annotate("segment",x=.42,xend=.57,y=.20,yend=.20,arrow=arrow(length=unit(1.7,"mm")),linewidth=.4,colour=C_GOLD) +
  annotate("segment",x=.80,xend=.72,y=.36,yend=.27,arrow=arrow(length=unit(1.7,"mm")),linewidth=.4,colour=C_GOLD) +
  annotate("text",x=.5,y=.06,label="Functional annotation is applied only after candidate selection.",size=2.2,colour=C_GRAY) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,0.86),clip="off") + labs(title="Cross-species evidence framework") + blank +
  theme(plot.title=element_text(size=8,face="bold",colour=C_TEXT,hjust=0))
fig1c <- tag(fig1c,"c")

# ---------------- Fig 2 ----------------
vst <- read_csv(file.path(ROOT,"06_meta_analysis/apple_main36/vst_biological_matrix.csv"), show_col_types=FALSE)
gene_col <- names(vst)[1]; mat <- as.matrix(vst[,-1]); rownames(mat) <- vst[[gene_col]]
pca <- prcomp(t(mat), center=TRUE, scale.=FALSE)
pc <- as.data.frame(pca$x[,1:2]); pc$sample <- rownames(pc)
pc$stage <- ifelse(grepl("_S1_",pc$sample),"Stage 1","Stage 2"); pc$rootstock <- ifelse(grepl("^M9_",pc$sample),"M9","MM106")
ve <- 100*(pca$sdev^2/sum(pca$sdev^2))
fig2a <- ggplot(pc,aes(PC1,PC2,shape=rootstock,colour=stage,fill=stage)) +
  geom_point(size=2.7, stroke=.75) +
  scale_colour_manual(values=c("Stage 1"=C_MINT,"Stage 2"="#ded4a3")) +
  scale_fill_manual(values=c("Stage 1"=C_MINT_LIGHT,"Stage 2"=C_YELLOW)) +
  scale_shape_manual(values=c("M9"=21,"MM106"=24)) +
  labs(x=sprintf("PC1 (%.1f%%)",ve[1]),y=sprintf("PC2 (%.1f%%)",ve[2]),colour="Stage",fill="Stage",shape="Rootstock") +
  guides(fill="none", colour=guide_legend(order=1,override.aes=list(shape=21,fill=c(C_MINT_LIGHT,C_YELLOW))), shape=guide_legend(order=2,override.aes=list(colour=C_GRAY,fill="white"))) +
  nature_theme + theme(legend.position="bottom")
fig2a <- tag(fig2a,"a")

de <- read_csv(file.path(ROOT,"06_meta_analysis/apple_main36/DE_summary.csv"),show_col_types=FALSE)
names(de)[1] <- "contrast"; de$contrast <- recode(de$contrast, stage1="Stage 1",stage2="Stage 2",interaction="Interaction")
fig2b <- ggplot(de,aes(x=reorder(contrast,FDR05),y=FDR05,fill=contrast)) + geom_col(width=.58) + coord_flip() +
  geom_text(aes(label=comma(FDR05)),hjust=-.12,size=2.35) + scale_fill_manual(values=c("Stage 1"=C_MINT,"Stage 2"=C_YELLOW,"Interaction"=C_PINK),guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,.15))) + labs(x=NULL,y="Genes (FDR < 0.05)") + nature_theme
fig2b <- tag(fig2b,"b")

# Fig 2c pathway concordance
k <- read_csv(file.path(ROOT,"06_meta_analysis/apple_main36/cross_stage_KEGG_concordance.csv"),show_col_types=FALSE) %>%
  filter(p.adjust_s1 < 0.1, p.adjust_s2 < 0.1) %>% mutate(class=ifelse(same_direction,"Same direction","Reversed"))
fig2c <- ggplot(k,aes(NES_s1,NES_s2,colour=class,fill=class)) +
  geom_hline(yintercept=0,colour=C_LIGHT,linewidth=.35) + geom_vline(xintercept=0,colour=C_LIGHT,linewidth=.35) +
  geom_abline(slope=1,intercept=0,linetype=2,colour=C_GRAY,linewidth=.35) +
  geom_point(shape=21,size=2.1,stroke=.62) +
  scale_colour_manual(values=c("Same direction"=C_MINT,"Reversed"="#e7bfd0")) +
  scale_fill_manual(values=c("Same direction"=C_MINT_LIGHT,"Reversed"="#f7dce8")) +
  annotate("text",x=-Inf,y=Inf,label="Conditional label-permutation tail fraction = 0.395",hjust=-.05,vjust=1.15,size=2.15,colour=C_GRAY) +
  labs(x="Stage 1 NES",y="Stage 2 NES") + nature_theme + theme(legend.position="bottom")
fig2c <- tag(fig2c,"c")

# Fig 2d external apple concordance
a <- read_csv(file.path(ROOT,"06_meta_analysis/apple2026/DE_new_shoot_PYGD076_vs_M9T337.csv"),show_col_types=FALSE) %>% select(GENEID,new_fc=log2FoldChange,new_padj=padj)
b <- read_csv(file.path(ROOT,"06_meta_analysis/apple2026/DE_fruiting_spur_PYGD076_vs_M9T337.csv"),show_col_types=FALSE) %>% select(GENEID,spur_fc=log2FoldChange,spur_padj=padj)
ext <- inner_join(a,b,by="GENEID") %>% filter(new_padj<.05 | spur_padj<.05) %>% mutate(class=ifelse(sign(new_fc)==sign(spur_fc),"Same direction","Discordant"))
pct_same <- mean(ext$class=="Same direction")*100
fig2d <- ggplot(ext,aes(new_fc,spur_fc,colour=class)) +
  geom_hline(yintercept=0,colour=C_LIGHT,linewidth=.35) + geom_vline(xintercept=0,colour=C_LIGHT,linewidth=.35) +
  geom_abline(slope=1,intercept=0,linetype=2,colour=C_GRAY,linewidth=.35) +
  geom_point(size=.85,alpha=.55) +
  scale_colour_manual(values=c("Same direction"=C_MINT,"Discordant"=C_PINK),guide="none") +
  annotate("label",x=-Inf,y=Inf,label=sprintf("%s genes significant in >=1 tissue\n%.1f%% same direction",comma(nrow(ext)),pct_same),hjust=-.03,vjust=1.08,size=2.15,linewidth=.2,fill="white",colour=C_TEXT) +
  labs(x="New-shoot leaves log2 fold change",y="Fruiting-spur leaves log2 fold change") + nature_theme
fig2d <- tag(fig2d,"d")


# Recompute from unmodified sample quantifications and persisted DEG/GSEA tables.
stopifnot(nrow(k)==37L, sum(k$class=="Reversed")==27L)
message("VALIDATED Fig2: pathways=",nrow(k), " reversed=",sum(k$class=="Reversed"),
        ", significant external genes=",nrow(ext), " pct same direction=",round(pct_same,1))
F2 <- (fig2a | fig2b) / (fig2c | fig2d)
save_both(fig2a, "Fig2a_PCA", 86, 66)
save_both(fig2b, "Fig2b_DEG_counts", 86, 66)
save_both(fig2c, "Fig2c_pathway_concordance", 86, 66)
save_both(fig2d, "Fig2d_external_apple_leaves", 86, 66)
save_both(F2, "Figure2_REVISED", 180, 142, COMB)
ggsave(file.path(COMB, "Figure2_REVISED.svg"), F2, width=180/25.4, height=142/25.4,
       units="in", device="svg", bg="white")
message("FIG2_RENDER_SUCCESS")
