suppressMessages(library(docopt))
library(plyr)
"Convert counts to RPKM/TPM for SE reads
Usage:
  counts_to_rpkm.R --counts=<counts> --gene_lengths=<gene_lengths> --outprefix=<outprefix> --inprefix=<inprefix> --gene_map=<gene_map>

  Options:
    --counts=<counts>                Counts files
    --inprefix=<inprefix>            Prefix for counts file
    --outprefix=<outprefix>          Outprefix
    --gene_lengths=<gene_lengths>    Gene lengths 
    --gene_map=<gene_map>            Filter for selecting only protein_coding genes

" -> doc
opts <- docopt(doc)
counts <- unlist(strsplit(opts[['counts']], split=','))
gene_lengths <- opts[['gene_lengths']]
outprefix <- opts[['outprefix']]
inprefix <- unlist(strsplit(opts[['inprefix']], split = ','))
gene_map <- opts[['gene_map']]

rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e6

}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6

}

countsList <- c()
gene_map <- read.table(gene_map, stringsAsFactors=F, header=F, row.names=1)
colnames(gene_map) <- c('gene_name', 'gene_type')
counts.df <- NULL
for (count in counts){
  count.df <- read.table(count, sep="\t", header=F, stringsAsFactors = F)
if (is.null(counts.df)){
    counts.df <- count.df
}
else{
  counts.df <- cbind(counts.df, count.df)
}
}
gene_lengths <- read.table(gene_lengths, stringsAsFactors=F, header=T)
#print(head(countsList[1]))
#counts.df <- cbind.fill(unlist(countsList))
#counts.df <- cbind(countsList[1], countsList[2])
#counts.df<- cbind(counts.WT, counts.3R, counts.50R)
counts.df<- counts.df[, c(1, seq(2, 2*length(counts),2))]
#df <- ldply(listOfDataFrames, data.frame)

names(counts.df)<- c("gene_name", inprefix)
#counts.df$gene_type <- gene_map[counts.df$gene_name,]$gene_type
#print(gene_map[row.names(counts.df),]$gene_type)i
lengths <- subset(gene_lengths, counts.df$gene_name %in% gene_lengths$Ensembl_gene_identifier)
lengths <- lengths[match(counts.df$gene_name, lengths$Ensembl_gene_identifier),]
counts.df$length <- lengths$length

write.table(counts.df, file=file.path(outprefix, paste('masterTable', 'counts', 'tsv', sep='.')), col.names = T)

for (col in inprefix){
  tpmcol <- tpm(counts.df[,col], counts.df$length)
  write.table(tpmcol, file=file.path(outprefix, paste(col, 'tpm', 'tsv', sep='.')), col.names = F)
}

for (col in inprefix){
  rpkmcol <- rpkm(counts.df[,col], counts.df$length)
  write.table(rpkmcol, file=file.path(outprefix, paste(col, 'rpkm', 'tsv', sep='.')), col.names = F)
}


