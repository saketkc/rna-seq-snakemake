library(EnsDb.Hsapiens.v86)
library(tximport)

dir <- ''

txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])

files <- files.path(dir, list.dirs(dir, recursive = F, full.names = F), 'counts',  'abundance.h5')
names(files) <- list.dirs(dir, recursive = F, full.names = F)
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene,  ignoreTxVersion = T)
write.table(txi.kallisto.tsv, 'merged_df.tsv', row.names=F)


