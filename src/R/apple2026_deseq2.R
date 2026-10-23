suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
})
root <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
outdir <- file.path(root,'06_meta_analysis','apple2026')
figdir <- file.path(root,'08_figures','apple2026')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)
meta <- read.delim(file.path(root,'01_metadata','apple2026_samples.tsv'), stringsAsFactors=FALSE)
meta$rootstock <- factor(meta$rootstock, levels=c('M9T337','PYGD076'))
meta$leaf_type <- factor(meta$leaf_type, levels=c('new_shoot','fruiting_spur'))
rownames(meta) <- meta$sample
files <- file.path(root,'04_quant','apple2026',meta$sample,'quant.sf')
names(files) <- meta$sample
stopifnot(all(file.exists(files)))
map_path <- file.path(root,"12_transcript_to_gene_mapping","apple_GCF_002114115.1_tx2gene_annot.csv")
if(file.exists(map_path)){
  map <- read.csv(map_path,stringsAsFactors=FALSE)
} else {
  ft <- read.delim(gzfile(file.path(root,'03_reference','apple_GDDH13','GCF_002114115.1_ASM211411v1_feature_table.txt.gz')), comment.char='', check.names=FALSE, stringsAsFactors=FALSE)
  colnames(ft)[1] <- sub('^# ','',colnames(ft)[1])
  map <- ft[ft$feature %in% c('mRNA','ncRNA') & nzchar(ft$product_accession) & nzchar(ft$GeneID), c('product_accession','GeneID','symbol','name')]
  map <- map[!duplicated(map$product_accession),]
  colnames(map)[1:2] <- c('TXNAME','GENEID')
}
stopifnot(all(c("TXNAME","GENEID","symbol","name") %in% names(map)))
tx <- tximport(files, type='salmon', tx2gene=map[,c('TXNAME','GENEID')], countsFromAbundance='lengthScaledTPM')
write.csv(data.frame(sample=colnames(tx$counts), library_size=colSums(tx$counts)), file.path(outdir,'library_sizes.csv'), row.names=FALSE)
dds <- DESeqDataSetFromTximport(tx, colData=meta, design=~ leaf_type + rootstock + leaf_type:rootstock)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
write.csv(as.data.frame(colData(dds)), file.path(outdir,'sample_design.csv'))
write.csv(as.data.frame(assay(vsd)), file.path(outdir,'vst_matrix.csv'))
res_new <- results(dds, name='rootstock_PYGD076_vs_M9T337', alpha=0.05)
int_name <- grep('leaf_typefruiting_spur.rootstockPYGD076', resultsNames(dds), value=TRUE)
stopifnot(length(int_name)==1)
res_int <- results(dds, name=int_name, alpha=0.05)
res_fruit <- results(dds, contrast=list(c('rootstock_PYGD076_vs_M9T337', int_name)), alpha=0.05)
annot <- unique(map[,c('GENEID','symbol','name')])
annot <- annot[!duplicated(annot$GENEID),]
export_res <- function(res, file){
  x <- as.data.frame(res); x$GENEID <- rownames(x)
  x <- merge(x, annot, by='GENEID', all.x=TRUE, sort=FALSE)
  x <- x[order(x$padj, na.last=TRUE),]
  write.csv(x, file, row.names=FALSE)
  x
}
x_new <- export_res(res_new, file.path(outdir,'DE_new_shoot_PYGD076_vs_M9T337.csv'))
x_fruit <- export_res(res_fruit, file.path(outdir,'DE_fruiting_spur_PYGD076_vs_M9T337.csv'))
x_int <- export_res(res_int, file.path(outdir,'DE_rootstock_by_tissue_interaction.csv'))
sig_count <- function(x){c(total=sum(!is.na(x$padj)), FDR05=sum(x$padj<0.05,na.rm=TRUE), FDR05_LFC1=sum(x$padj<0.05 & abs(x$log2FoldChange)>=1,na.rm=TRUE))}
sumtab <- rbind(new_shoot=sig_count(res_new), fruiting_spur=sig_count(res_fruit), interaction=sig_count(res_int))
write.csv(sumtab, file.path(outdir,'DE_summary.csv'))
pca <- plotPCA(vsd, intgroup=c('rootstock','leaf_type'), returnData=TRUE)
percentVar <- round(100 * attr(pca,'percentVar'))
p <- ggplot(pca, aes(PC1,PC2,color=rootstock,shape=leaf_type,label=name)) + geom_point(size=3) + geom_text(vjust=-0.8, size=2.5) + xlab(paste0('PC1: ',percentVar[1],'%')) + ylab(paste0('PC2: ',percentVar[2],'%')) + theme_bw()
ggsave(file.path(figdir,'Fig_PCA_apple2026.png'), p, width=7, height=5, dpi=300)
volcano <- function(x,title,file){
  x$neglog10 <- -log10(pmax(x$padj,1e-300)); x$sig <- !is.na(x$padj) & x$padj<0.05 & abs(x$log2FoldChange)>=1
  q <- ggplot(x,aes(log2FoldChange,neglog10)) + geom_point(aes(alpha=sig),size=1) + scale_alpha_manual(values=c('FALSE'=0.25,'TRUE'=0.9),guide='none') + geom_vline(xintercept=c(-1,1),linetype=2) + geom_hline(yintercept=-log10(0.05),linetype=2) + labs(title=title,x='log2 fold change',y='-log10 adjusted P') + theme_bw()
  ggsave(file,q,width=6,height=5,dpi=300)
}
volcano(x_new,'New shoot: PYG D076 vs M9T337',file.path(figdir,'Volcano_new_shoot.png'))
volcano(x_fruit,'Fruiting spur: PYG D076 vs M9T337',file.path(figdir,'Volcano_fruiting_spur.png'))
cat('RESULTS_NAMES\n'); print(resultsNames(dds)); cat('SUMMARY\n'); print(sumtab)
saveRDS(dds,file.path(outdir,'dds.rds')); saveRDS(vsd,file.path(outdir,'vsd.rds'))
