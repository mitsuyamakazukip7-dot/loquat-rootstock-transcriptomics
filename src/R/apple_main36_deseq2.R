suppressPackageStartupMessages({library(tximport); library(DESeq2); library(ggplot2); library(clusterProfiler)})
root <- path.expand(Sys.getenv("LOQUAT_PROJECT_ROOT",unset="~/Desktop/LOQUAT_DWARFING_PUBLIC"))
outdir <- file.path(root,'06_meta_analysis','apple_main36'); figdir <- file.path(root,'08_figures','apple_main36')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE); dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
meta <- read.delim(file.path(root,'01_metadata','apple_main36_samples.tsv'),stringsAsFactors=FALSE)
qfiles <- file.path(root,'04_quant','apple_main36',meta$sample,'quant.sf'); names(qfiles) <- meta$sample
stopifnot(length(qfiles)==36, all(file.exists(qfiles)))
map_path <- file.path(root,"12_transcript_to_gene_mapping","apple_GCF_002114115.1_tx2gene_annot.csv")
if(file.exists(map_path)){
  map <- read.csv(map_path,stringsAsFactors=FALSE)
} else {
  ft <- read.delim(gzfile(file.path(root,'03_reference','apple_GDDH13','GCF_002114115.1_ASM211411v1_feature_table.txt.gz')),comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
  colnames(ft)[1] <- sub('^# ','',colnames(ft)[1])
  map <- ft[ft$feature %in% c('mRNA','ncRNA') & nzchar(ft$product_accession) & nzchar(ft$GeneID),c('product_accession','GeneID','symbol','name')]
  map <- map[!duplicated(map$product_accession),]; colnames(map)[1:2] <- c('TXNAME','GENEID')
}
stopifnot(all(c("TXNAME","GENEID","symbol","name") %in% names(map)))
tx <- tximport(qfiles,type='salmon',tx2gene=map[,c('TXNAME','GENEID')],countsFromAbundance='lengthScaledTPM')
meta$bio_rep <- sub('-.*$','',meta$rep); meta$tech_rep <- sub('^.*-','',meta$rep)
meta$bio_id <- paste(meta$rootstock,paste0('S',meta$stage),paste0('B',meta$bio_rep),sep='_')
# conservative pseudo-replication control: sum 3 technical/library replicates within biological replicate
bi <- unique(meta[,c('bio_id','rootstock','stage','bio_rep')]); rownames(bi) <- bi$bio_id
agg <- sapply(bi$bio_id,function(id) rowSums(tx$counts[,meta$bio_id==id,drop=FALSE]))
colnames(agg) <- bi$bio_id
write.csv(data.frame(sample=meta$sample,bio_id=meta$bio_id,rootstock=meta$rootstock,stage=meta$stage,bio_rep=meta$bio_rep,tech_rep=meta$tech_rep),file.path(outdir,'technical_to_biological_map.csv'),row.names=FALSE)
write.csv(data.frame(bio_id=colnames(agg),library_size=colSums(agg)),file.path(outdir,'aggregated_library_sizes.csv'),row.names=FALSE)
bi$rootstock <- factor(bi$rootstock,levels=c('MM106','M9')); bi$stage <- factor(bi$stage,levels=c('1','2'))
dds <- DESeqDataSetFromMatrix(round(agg),colData=bi,design=~stage+rootstock+stage:rootstock)
dds <- dds[rowSums(counts(dds)>=10)>=3,]; dds <- DESeq(dds)
vsd <- vst(dds,blind=FALSE); write.csv(as.data.frame(assay(vsd)),file.path(outdir,'vst_biological_matrix.csv'))
res_s1 <- results(dds,name='rootstock_M9_vs_MM106',alpha=0.05)
intn <- grep('stage2.rootstockM9',resultsNames(dds),value=TRUE); stopifnot(length(intn)==1)
res_int <- results(dds,name=intn,alpha=0.05)
res_s2 <- results(dds,contrast=list(c('rootstock_M9_vs_MM106',intn)),alpha=0.05)
ann <- unique(map[,c('GENEID','symbol','name')]); ann <- ann[!duplicated(ann$GENEID),]
expres <- function(res,fn){x<-as.data.frame(res);x$GENEID<-rownames(x);x<-merge(x,ann,by='GENEID',all.x=TRUE,sort=FALSE);x<-x[order(x$padj,na.last=TRUE),];write.csv(x,file.path(outdir,fn),row.names=FALSE);x}
x1<-expres(res_s1,'DE_stage1_M9_vs_MM106.csv'); x2<-expres(res_s2,'DE_stage2_M9_vs_MM106.csv'); xi<-expres(res_int,'DE_stage_interaction.csv')
sc<-function(r)c(tested=sum(!is.na(r$padj)),FDR05=sum(r$padj<.05,na.rm=TRUE),FDR05_LFC1=sum(r$padj<.05&abs(r$log2FoldChange)>=1,na.rm=TRUE))
write.csv(rbind(stage1=sc(res_s1),stage2=sc(res_s2),interaction=sc(res_int)),file.path(outdir,'DE_summary.csv'))
pca<-plotPCA(vsd,intgroup=c('rootstock','stage'),returnData=TRUE);pv<-round(100*attr(pca,'percentVar'))
p<-ggplot(pca,aes(PC1,PC2,color=rootstock,shape=stage,label=name))+geom_point(size=3)+geom_text(vjust=-.8,size=2.5)+xlab(paste0('PC1: ',pv[1],'%'))+ylab(paste0('PC2: ',pv[2],'%'))+theme_bw();ggsave(file.path(figdir,'PCA_biological_replicates.png'),p,width=7,height=5,dpi=300)
# technical/library replicate QC before aggregation
rawdds <- DESeqDataSetFromMatrix(round(tx$counts),colData=meta,design=~1); rawdds<-rawdds[rowSums(counts(rawdds)>=10)>=3,]
rawv <- vst(rawdds,blind=TRUE); cm <- cor(assay(rawv),method='spearman')
pairs <- which(upper.tri(cm),arr.ind=TRUE); pq <- data.frame(s1=colnames(cm)[pairs[,1]],s2=colnames(cm)[pairs[,2]],rho=cm[pairs])
rownames(meta) <- meta$sample; pq$bio1 <- meta[pq$s1,'bio_id']; pq$bio2 <- meta[pq$s2,'bio_id']; pq$root1 <- meta[pq$s1,'rootstock']; pq$root2 <- meta[pq$s2,'rootstock']; pq$stage1 <- meta[pq$s1,'stage']; pq$stage2 <- meta[pq$s2,'stage']
pq$class <- ifelse(pq$bio1==pq$bio2,'within_biological_group',ifelse(pq$root1==pq$root2 & pq$stage1==pq$stage2,'between_bio_same_condition','other'))
write.csv(pq,file.path(outdir,'technical_rep_pairwise_correlations.csv'),row.names=FALSE)
qcmed<-aggregate(rho~class,pq,median); write.csv(qcmed,file.path(outdir,'technical_rep_correlation_summary.csv'),row.names=FALSE)
# cross-stage effect concordance
m<-merge(x1[,c('GENEID','log2FoldChange','padj','symbol','name')],x2[,c('GENEID','log2FoldChange','padj')],by='GENEID',suffixes=c('_s1','_s2'))
ok<-is.finite(m$log2FoldChange_s1)&is.finite(m$log2FoldChange_s2);m<-m[ok,]
either<-(!is.na(m$padj_s1)&m$padj_s1<.05)|(!is.na(m$padj_s2)&m$padj_s2<.05);both<-(!is.na(m$padj_s1)&m$padj_s1<.05)&(!is.na(m$padj_s2)&m$padj_s2<.05)
metrics<-data.frame(metric=c('n','spearman','pearson','same_direction_all','either_sig_n','same_direction_either','both_sig_same_direction_n'),value=c(nrow(m),cor(m$log2FoldChange_s1,m$log2FoldChange_s2,method='spearman'),cor(m$log2FoldChange_s1,m$log2FoldChange_s2),mean(sign(m$log2FoldChange_s1)==sign(m$log2FoldChange_s2)),sum(either),mean(sign(m$log2FoldChange_s1[either])==sign(m$log2FoldChange_s2[either])),sum(both & sign(m$log2FoldChange_s1)==sign(m$log2FoldChange_s2))))
write.csv(metrics,file.path(outdir,'cross_stage_concordance.csv'),row.names=FALSE); write.csv(m[both & sign(m$log2FoldChange_s1)==sign(m$log2FoldChange_s2),],file.path(outdir,'shared_DE_same_direction.csv'),row.names=FALSE)
# Canonical rerun uses the fixed full-gene KEGG map; the historical gseKEGG
# output (37 pathways, 27 reversed) remains archived as original provenance.
# Frozen rerun is a sensitivity calculation and may return 38/28.
kegg_path <- file.path(root,"07_KEGG_GSEA_full_gene_universe","KEGG_MDM_FULL_TERM2GENE_20261022.csv")
stopifnot(file.exists(kegg_path))
frozen_kegg<-read.csv(kegg_path,stringsAsFactors=FALSE)
stopifnot(all(c("pathway","gene") %in% names(frozen_kegg)))
gsea_kegg <- function(x,prefix){
  y<-x[is.finite(x$stat)&!is.na(x$GENEID),]; y$GENEID<-as.character(y$GENEID); y<-y[!duplicated(y$GENEID),]
  gl<-y$stat; names(gl)<-y$GENEID; gl<-sort(gl,decreasing=TRUE)
  set.seed(if(prefix=="stage1")20261022 else 20261023)
  g<-clusterProfiler::GSEA(geneList=gl,TERM2GENE=frozen_kegg[,c("pathway","gene")],
      minGSSize=10,maxGSSize=500,pvalueCutoff=1,verbose=FALSE,seed=TRUE,by="fgsea")
  z<-as.data.frame(g); write.csv(z,file.path(outdir,paste0(prefix,'_GSEA_FULL_FROZEN.csv')),row.names=FALSE); z
}
g1<-gsea_kegg(x1,'stage1'); g2<-gsea_kegg(x2,'stage2')
kg<-merge(g1[,c('ID','Description','NES','p.adjust')],g2[,c('ID','NES','p.adjust')],by='ID',suffixes=c('_s1','_s2'))
kg$same_direction<-sign(kg$NES_s1)==sign(kg$NES_s2);kg$maxFDR<-pmax(kg$p.adjust_s1,kg$p.adjust_s2);kg<-kg[order(kg$maxFDR),]
write.csv(kg,file.path(outdir,'cross_stage_KEGG_concordance.csv'),row.names=FALSE)
saveRDS(dds,file.path(outdir,'dds_biological_aggregated.rds'));saveRDS(vsd,file.path(outdir,'vsd_biological_aggregated.rds'))
cat('RESULTS_NAMES\n');print(resultsNames(dds));cat('QC_TECH_MEDIAN\n');print(qcmed);cat('DE_SUMMARY\n');print(rbind(stage1=sc(res_s1),stage2=sc(res_s2),interaction=sc(res_int)))
cat('CROSS_STAGE\n');print(metrics);cat('KEGG_SHARED_FDR10_SAME_DIR=',sum(kg$same_direction&kg$p.adjust_s1<.1&kg$p.adjust_s2<.1,na.rm=TRUE),'\n')
