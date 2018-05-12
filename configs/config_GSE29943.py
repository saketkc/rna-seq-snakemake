GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/home/cmb-panasas2/skchoudh/public_HuR/GSE29943/rna_seq'
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/public_HuR/GSE29943'
SAMPLES=['SRR309282', 'SRR309283',
         'SRR309284']#, 'SRR309288',
         #'SRR309289', 'SRR309290']

GENOME_BUILD = 'hg38'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/rna-seq-snakemake/scripts'
STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.annotation.gtf'

GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
GENE_NAME_MAP = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD + '_gene_names_stripped.tsv'



