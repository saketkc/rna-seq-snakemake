include:
    'config_cluster.py'

workdir: OUT_DIR

from itertools import chain
from os.path import join

rule all:
    input:
        STAR_INDEX,
        expand('qc/{sample}_R1_001_fastqc.html', sample=SAMPLES),
        expand('qc/{sample}_R2_001_fastqc.html', sample=SAMPLES),
        expand('mapped/counts_strict/star/{sample}.counts.tsv', sample=SAMPLES)


rule create_index:
    input:
        GENOME_FASTA,
        GTF
    output: STAR_INDEX
    threads: 16
    shell:
        r'''mkdir -p {output} && STAR --runThreadN 16\
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input[0]}\
            --sjdbGTFfile {input[1]}'''


rule perform_qc:
    input:
        R1=RAWDATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAWDATA_DIR+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir = 'qc'
    output:
       'qc/{sample}_R1_001_fastqc.html',
       'qc/{sample}_R1_001_fastqc.zip',
       'qc/{sample}_R2_001_fastqc.html',
       'qc/{sample}_R2_001_fastqc.zip',
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input.R1} {input.R2}
        '''

rule perfom_trimming:
    input:
        R1=RAWDATA_DIR+'/{sample}_R1_001.fastq.gz',
        R2=RAWDATA_DIR+'/{sample}_R2_001.fastq.gz'
    params:
        out_dir='preprocessed/{sample}',
        phred_cutoff=5
    output:
        'preprocessed/{sample}/{sample}_R1_001_val_1.fq.gz',
        'preprocessed/{sample}/{sample}_R2_001_val_2.fq.gz',
    shell:
        r'''
            trim_galore --paired -o {params.out_dir} -q {params.phred_cutoff} {input.R1} {input.R2}
        '''

"""
rule map_star:
    input:
        R1='preprocessed/{sample}/{sample}_R1_001_val_1.fq.gz',
        R2='preprocessed/{sample}/{sample}_R2_001_val_2.fq.gz',
        index=STAR_INDEX
    output: 'mapped/bams/star/{sample}.bam'
    params:
        prefix = 'mapped/bams/star/{sample}',
        unmapped = 'unmapped/fastq/star/{sample}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STAR --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --outFilterMatchNmin 50\
             --outFilterMismatchNmax 100\
             --readFilesCommand zcat\
             --outReadsUnmapped {params.unmapped} && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
        '''
"""

rule map_starlong:
    input:
        R1='preprocessed/{sample}/{sample}_R1_001_val_1.fq.gz',
        R2='preprocessed/{sample}/{sample}_R2_001_val_2.fq.gz',
        index=STAR_INDEX
    output: 'mapped/bams/star/{sample}.bam'
    params:
        prefix = 'mapped/bams/star/{sample}',
        unmapped = 'unmapped/fastq/star/{sample}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STARlong --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1} {input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --readFilesCommand zcat\
             --outReadsUnmapped {params.unmapped}\
            --outFilterMultimapScoreRange 20\
            --outFilterScoreMinOverLread 0\
            --outFilterMatchNminOverLread 0.66\
            --outFilterMismatchNmax 100\
            --winAnchorMultimapNmax 200\
            --seedPerReadNmax 100000 --seedPerWindowNmax 100\
            && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} &&\
            mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out\
            {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
        '''
rule sort_by_name:
    input: 'mapped/bams/star/{sample}.bam'
    output: 'mapped/bams/star/{sample}.sortedByName.bam'
    shell:
        r'''
            samtools sort -on {input} -T /tmp/ -o {output}
        '''
rule count:
    input: 'mapped/bams/star/{sample}.sortedByName.bam'
    params:
        annotation=GTF,
        phred_cutoff=5
    output: 'mapped/counts_strict/star/{sample}.counts.tsv'
    shell:
        r'''
        source activate clipseq2 && htseq-count --order=name --format=bam --mode=intersection-strict --stranded=no --minaqual={params.phred_cutoff} --type=exon --idattr=gene_id {input} {params.annotation} > {output}
        '''
