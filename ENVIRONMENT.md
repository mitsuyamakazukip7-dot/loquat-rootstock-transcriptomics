# Tested environment, October 2026

Local research Mac: R 4.5.1; DESeq2 1.50.2, tximport 1.38.2, clusterProfiler 4.18.4, fgsea 1.36.2, limma 3.66.0, KEGGREST 1.50.0.
ggplot2 4.0.3; patchwork 1.3.2; svglite 2.2.2; plus dplyr, tidyr, readr, scales.
Python 3.14.7; NumPy 2.5.3, Matplotlib 3.11.2, pandas 3.0.6 from the local project virtual environment; openpyxl required only for optional XLSX builder.

Public-read reprocessing additionally needs SRA Toolkit, Salmon 2.7.0, and reference Salmon indices; raw protein RBH needs DIAMOND and original FASTA.

Source-specific source tree paths/indices must be configured for upstream reanalysis. The frozen pathway GSEA/permutation has its own small standalone inputs. Mac success does not establish clean-machine full-FASTQ E2E portability.