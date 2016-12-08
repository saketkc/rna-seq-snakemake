GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/staging/as/skchoudh/rna/GSE62129'
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/data/public_HuR/GSE62129'
GENOME_BUILD = 'mm10'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/rna-seq-snakemake/scripts'
STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.annotation.gtf'

GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
GENE_NAME_MAP = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD + '_gene_names_stripped.tsv'
DESIGN_FILE = RAWDATA_DIR + '/design.txt'


