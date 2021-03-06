GENOMES_DIR='/home/cmb-panasas2/skchoudh/genomes'
OUT_DIR = '/staging/as/skchoudh/rna/HuR_results/mouse/rna_seq_star_mm10_annotated'
#OUT_DIR = '/home/cmb-06/as/skchoudh/rna/HuR_results/mouse/rna_seq_star_mm10_annotated/'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/clip_seq_pipeline/scripts'
RAWDATA_DIR ='/home/cmb-06/as/skchoudh/dna/HuR_Mouse_Human_liver/rna-seq/Penalva_L_08182016/mouse'
GENOME_BUILD = 'mm10'
GENOME_FASTA = GENOMES_DIR + '/' + GENOME_BUILD + '/fasta/'+ GENOME_BUILD+ '.fa'

STAR_INDEX = GENOMES_DIR + '/' + GENOME_BUILD + '/star_annotated'
GTF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.annotation.gtf'
GENE_NAMES = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + GENOME_BUILD+'_gene_names_stripped.tsv'
GENE_LENGTHS = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.coding_lengths.tsv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
DESIGN_FILE = '/home/cmb-06/as/skchoudh/dna/HuR_Mouse_Human_liver/rna-seq/Penalva_L_08182016/exp-design/mouse_design.txt'
INTERGENIC_BED = '/staging/as/skchoudh/rna/HuR_results/conservation/annotated_beds/union/mmaster_confident.annotated.intergenic.bed'
ANNOT_DIR = '/home/cmb-panasas2/skchoudh/github_projects/gencode_regions/data/GRCm38/vM11/'
GENE_BED = ANNOT_DIR + 'genes.bed'

HTSEQ_STRANDED = 'no'
FEATURECOUNTS_S = ''

GENE_BED = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.annotation.genePred.ssv'  #+ GENOME_BUILD+'_gene_lengths.tsv'
GENE_BED_TSV = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.annotation.genePred'  #+ GENOME_BUILD+'_gene_lengths.tsv'
DEXSeq_GFF = GENOMES_DIR + '/' + GENOME_BUILD + '/annotation/' + 'gencode.vM11.annotation.DEXSeq.gff'
