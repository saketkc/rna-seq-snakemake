#!/usr/bin/env Rscript

suppressMessages(library(docopt))
suppressMessages(library(readr))
#suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
#suppressMessages(library(limma))
suppressMessages(library('RColorBrewer'))
suppressMessages(library('gplots'))

suppressMessages(library('BiocParallel'))
register(MulticoreParam(4))

"Perform differetioan expressiona analysis
Usage:
  do_DE_analysis.R --basedir=<basedir> --design_file=<design_file> --gene_annotations=<gene_annotations> --outprefix=<outprefix> --inprefix=<inprefix>

  Options:
    --basedir=<basedir>                      Absolute path to base directory
    --design_file=<design_file>              Absolute path to design matrix file
    --inprefix=<inprefix>                    Prefix for counts file
    --outprefix=<outprefix>                  Outprefix
    --gene_annotations=<gene_annotations>    Annotation file

" -> doc

opts <- docopt(doc)

base_dir <- opts[['basedir']]   #'/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/mapped/counts_strict/star'
design_file <- opts[['design_file']] #'/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/mapped/mouse_design.txt'
outprefix <- opts[['outprefix']]  #'/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/DE_analysis/mouse_de_analysis'
gene_annotations <- opts[['gene_annotations']] #'/media/dna/genomes/mm10/annotation/mm10_gene_names.tsv'
gene_annotations <- read.table(gene_annotations, row.names = 1, col.names = c('id', 'name', 'type'))
inprefix <- opts[['inprefix']]

#base_dir <- '/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/mapped/counts_strict/star'
#design_file <- '/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/mapped/mouse_design.txt'
#outprefix <- '/media/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/DE_analysis/mouse_de_analysis'
#gene_annotations <- '/media/dna/genomes/mm10/annotation/mm10_gene_names.tsv'
#inprefix <- 'counts'

#base_dir <- '/media/rna/HuR_results/human/rna_seq_star_hg38_annotated/mapped/counts_strict/star'
#design_file <- '/media/rna/HuR_results/human/rna_seq_star_hg38_annotated/mapped/human_design_chucked.txt'
#outprefix <- '/media/rna/HuR_results/human/rna_seq_star_hg38_annotated/DE_analysis/human_de_analysis'
#gene_annotations <- '/media/dna/genomes/hg38/annotation/hg38_gene_names.tsv'
#inprefix <- 'counts'

#gene_annotations <- read.table(gene_annotations, row.names = 1, col.names = c('id', 'name', 'type'))

design.info <- read.csv(design_file, header = TRUE, stringsAsFactors=FALSE)
sample_id <- design.info$sample
files <- paste(sample_id, inprefix, 'tsv', sep='.')
names(files) <- sample_id
condition <- design.info$condition

sampleTable <- data.frame(sampleName=sample_id, fileName=files, condition = condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=base_dir, design=~condition)
## We dont need version numbers
rownames(ddsHTSeq) <- gsub('\\.[0-9]+', '', rownames(ddsHTSeq))
## Filter genes with atleast 2 count
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1,  ]

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','knockdown'))

## Analyze with DESeq2
dds <- DESeq(ddsHTSeq)
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), paste(outprefix, 'DESeq2', 'all', 'tsv', sep='.'))
resSig <- subset(resOrdered, padj < 0.01)
#row.names(resSig) <- gene_annotations[row.names(resSig),]$name
write.table(as.data.frame(resSig), paste(outprefix, 'DESeq2', 'sig', 'tsv', sep='.'))

row.names(resSig) <- gene_annotations[row.names(resSig),]$name
write.table(as.data.frame(resSig), paste(outprefix, 'DESeq2', 'sig', 'names', 'tsv', sep='.'))
