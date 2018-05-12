GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/staging/as/skchoudh/rna/public_HuR/JBC_UREN_2011/rna-seq'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'
RAWDATA_DIR ='/staging/as/skchoudh/dna/public_HuR/JBC_UREN_2011/rna-seq'
GENOME_BUILD = 'hg38'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'
STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.annotation.gtf'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
DESIGN_FILE = '/staging/as/skchoudh/dna/Public_HuR/rna-seq/design.txt'
ANNOT_DIR = '/home/cmb-panasas2/skchoudh/github_projects/gencode_regions/data/GRCh37/v25/'
#GENE_BED = ANNOT_DIR + 'genes.bed'

GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.24.genes.bed'  #+ GENOME_BUILD+'_gene_lengths.tsv'
