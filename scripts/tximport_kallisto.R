library(EnsDb.Hsapiens.v86)
library(tximport)
library(DESeq2)


dir <- ''

txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])

files <- files.path(dir, list.dirs(dir, recursive = F, full.names = F), 'counts',  'abundance.h5')
names(files) <- list.dirs(dir, recursive = F, full.names = F)
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene,  ignoreTxVersion = T)

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

rld <- rlog(dds, blind=FALSE)

write.table(rld, 'merged_df.tsv', row.names=F)


