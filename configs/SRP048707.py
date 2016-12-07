####################################################Input##################################################

## Specify barcode length
BARCODE_LENGTH = 9

## Specify genome build and species of species to which the fastqs belong
GENOME_BUILD = 'mm10'
SPECIES = 'mouse'

## Specify genome build of other species for liftover
GENOME_BUILD_OTHER = 'hg38'

## Parent location of raw .fq files
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/data/public_HuR/GSE62129'

## Absolute path to the scripts directory
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'

## Absolute path to root directory where results will be created
ANALYSIS_DIR = '/staging/as/skchoudh/rna/GSE62129'
###########################################################################################################



##############################################GENOME specific##############################################
GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'

GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'fasta' + '/' + GENOME_BUILD + '.fa'

GENOME_SIZES = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'fasta' + '/' + GENOME_BUILD + '.sizes'

STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'star_annotated'

UTR5_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.5UTRs.bed'
UTR3_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.3UTRs.bed'
CDS_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.CDS.bed'
INTRON_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.introns.bed'
MIRNA_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'miRNA' + '/'+ 'GrCM38.miRNA.bed'
GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.genes.bed'
LIFTOVER_CHAIN = GENOMES_DIR + '/' + GENOME_BUILD + 'liftover' + '/' + 'mm10ToHg38.over.chain'

LINCRNA_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'lncRNA' + '/'+ 'gencode.vM11.long_noncoding_RNAs.named.bed'
GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + 'gencode.vM11.genes.bed'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/' + 'annotation' + '/' + GENOME_BUILD + '_gene_names.tsv'
###########################################################################################################


################################################OTHER GENOME################################################

UTR5_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.v25.5UTRs.bed'

UTR3_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.v25.3UTRs.bed'

CDS_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.v25.CDS.bed'

INTRON_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.v25.introns.bed'

MIRNA_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'miRNA' + '/' + 'GrCh38.miRNA.bed'

GENE_BED_OTHER = GENOMES_DIR + '/' + GENOME_BUILD_OTHER + '/' + 'annotation' + '/' + 'gencode.v25.genes.bed'

LIFTOVER_CHAIN_OTHER = GENOMES_DIR + '/' + GENOME_BUILD + 'liftover' + '/' + 'hg38ToMm10.over.chain'
###########################################################################################################


################################################DO NOT EDIT#################################################
LIFT_PREFIX = GENOME_BUILD + 'To' + GENOME_BUILD_OTHER
############################################################################################################
