CDNA = '/home/cmb-panasas2/skchoudh/genomes/hg19/kallisto/hg19'
GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/home/cmb-panasas2/skchoudh/HuR_results/human/rna_seq_star_hg38_annotated'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/data/HuR_Mouse_Human_liver/rna-seq/Penalva_L_08182016'
SAMPLES=['HepG2_CTRL1_S31_L004', 'HepG2_CTRL_2_S33_L004',
         'HepG2_CTRL_7_S35_L004', 'HepG2_HuR_KD_1_S32_L004',
         'HepG2_HuR_KD_2_S34_L004', 'HepG2_HuR_KD_7_S36_L004']

GENOME_BUILD = 'hg38'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'

STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.annotation.gtf'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
GENE_NAME_MAP = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD + '_gene_names_stripped.tsv'



