GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/staging/as/skchoudh/rna/HuR_results/human/rna_seq_Penalva_L_01182017'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'
RAWDATA_DIR ='/staging/as/skchoudh/dna/lai_data/Penalva_L_01182017'
SAMPLES=['HepG2_CTRL1_S31_L004', 'HepG2_CTRL_2_S33_L004',
         'HepG2_CTRL_7_S35_L004', 'HepG2_HuR_KD_1_S32_L004',
         'HepG2_HuR_KD_2_S34_L004', 'HepG2_HuR_KD_7_S36_L004']

GENOME_BUILD = 'hg38'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'

STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.annotation.gtf'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.v25.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'




