#' Convert counts to transcripts per million (TPM).
#' 
#' Convert a numeric matrix of features (rows) and conditions (columns) with
#' raw feature counts to transcripts per million.
#' 
#'    Lior Pachter. Models for transcript quantification from RNA-Seq.
#'    arXiv:1104.3889v2 
#'    
#'    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#'    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#'    doi:10.1007/s12064-012-0162-3
#'    
#' @param counts A numeric matrix of raw feature counts i.e.
#'  fragments assigned to each gene.
#' @param featureLength A numeric vector with feature lengths.
#' @param meanFragmentLength A numeric vector with mean fragment lengths.
#' @return tpm A numeric matrix normalized by library size and feature length.
#' "
suppressMessages(library(docopt))
library(plyr)
cbind.fill <- function(...) {                                                                                                                                                       
  transpoted <- lapply(list(...),t)                                                                                                                                                 
  transpoted_dataframe <- lapply(transpoted, as.data.frame)                                                                                                                         
  return (data.frame(t(rbind.fill(transpoted_dataframe))))                                                                                                                          
}
"Convert counts to TPM
Usage:
  counts_to_tpm.R --counts=<counts>... --insert_sizes=<insert_sizes>... --gene_lengths=<gene_lengths> --outprefix=<outprefix> --inprefix=<inprefix> --gene_map=<gene_map>

  Options:
    --counts=<counts>                Counts files
    --inprefix=<inprefix>            Prefix for counts file
    --outprefix=<outprefix>          Outprefix
    --gene_lengths=<gene_lengths>    Gene lengths 
    --insert_sizes=<insert_sizes>    File with insert sizes
    --gene_map=<gene_map>            Filter for selecting only protein_coding genes

" -> doc
opts <- docopt(doc)

counts <- unlist(strsplit(opts[['counts']], split=','))
insert_sizes <- unlist(strsplit(opts[['insert_sizes']], split=','))
gene_lengths <- opts[['gene_lengths']]
outprefix <- opts[['outprefix']]
inprefix <- unlist(strsplit(opts[['inprefix']], split = ','))
gene_map <- opts[['gene_map']]

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(meanFragmentLength) == ncol(counts))
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)[idx]
  return(tpm)
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
#print(gene_map[row.names(counts.df),]$gene_type)

#protein_coding.counts <- subset(counts.df, counts.df$gene_type == 'protein_coding')

protein_coding.counts <- subset(counts.df, counts.df$gene_name %in% gene_lengths$Ensembl_gene_identifier)
protein_coding.counts <- protein_coding.counts[match(gene_lengths$Ensembl_gene_identifier, protein_coding.counts$gene_name),]
#counts_df <- as.data.frame(counts.df)
insert_sizes <- as.numeric(unname(sapply(insert_sizes, readLines)))

rawcounts <- as.matrix(protein_coding.counts[, seq(2,length(counts)+1) ])

rownames(rawcounts)<- protein_coding.counts$gene_name
TPM <- counts_to_tpm(rawcounts, gene_lengths$length, insert_sizes)

for (col in inprefix){
  tpm <- TPM[, col]
  write.table(tpm, file=file.path(outprefix, paste(col, 'tpm', 'tsv', sep='.')), col.names = F)
}
