from itertools import chain, combinations
from os.path import join
import glob
import re
from collections import defaultdict
import os
from os.path import join
import glob
import numpy as np
import pandas as pd
import os
import errno

def mkdir_p(path):
  """Python version mkdir -p

  Parameters
  ----------

  path : str
  """
  if path:
    try:
      os.makedirs(path) 
    except OSError as exc:  # Python >2.5
      if exc.errno == errno.EEXIST and os.path.isdir(path):
        pass
      else:
        raise

include:
    config['config_path']

workdir: OUT_DIR


mkdir_p(join(OUT_DIR, 'slurm-logs'))
ALL_SRA_FILES = glob.glob('{}/**/*.sra'.format(RAWDATA_DIR), recursive=True)
SRX_ID = defaultdict(list)

for sample in ALL_SRA_FILES:
    srx, srr = sample.replace('{}'.format(RAWDATA_DIR),'').lstrip('/').rstrip('/').split('/')
    SRX_ID[srx].append(srr.replace('.sra', ''))

SRR_ID = list(SRX_ID.values())
SRX_SAMPLES = list(SRX_ID.keys())

ALL_SRR = [item for sublist in SRR_ID for item in sublist]


def merge_fastq_input(wildcards):
    return ['sratofastq/{}.fastq.gz'.format(srr) for srr in SRX_ID[wildcards.sample]]


def sra_to_fastq_input(wildcards):
    srr_id = wildcards.sample
    for srx_id in list(SRX_ID.keys()):
        value = SRX_ID[srx_id]
        if srr_id in list(value):
            return ancient(str(join(RAWDATA_DIR, srx_id, srr_id+'.sra')))
    print("WRONG encodeterend: {}".format(srr_id))

def get_wrapper(wrapper_name):
    path = os.path.dirname(os.path.abspath(os.path.realpath(workflow.snakefile)))
    return 'file://' + os.path.join(path,
                                    'wrappers',
                                    wrapper_name + '.py')


rule all:
    input:
        CDNA_IDX,
        CDS_IDX,
        expand('counts_cdna/{sample}/abundance.tsv', sample=SRX_SAMPLES),
        expand('counts_cds/{sample}/abundance.tsv', sample=SRX_SAMPLES)


rule sra_to_fastq:
    input: sra_to_fastq_input
    benchmark: 'benchmarks/sra_to_fastq/{sample}.txt'
    output: 'sratofastq/{sample}.fastq.gz'
    params:
        prefix_sra='sratofastq/{sample}.sra.fastq',
        prefix='sratofastq/{sample}.fastq',
        temp_dir='/tmp/{sample}_sratofastq',
    shell:
        r'''
        fastq-dump --split-3 \
        -O sratofastq {input} \
        && gzip {params.prefix}
        '''
"""

rule sra_to_fastq:
    input: sra_to_fastq_input
    benchmark: 'benchmarks/sra_to_fastq/{sample}.txt'
    output: 'sratofastq/{sample}.fastq.gz'
    params:
        prefix_sra='sratofastq/{sample}.sra.fastq',
        prefix='sratofastq/{sample}.fastq',
        temp_dir='/tmp/{sample}_sratofastq',
    threads: 16
    shell:
        r'''
        parallel-fastq-dump --threads 16 --outdir sratofastq/ --gzip -s {input}
        '''
"""

rule merge_fastq:
    input:
        all_fastq = expand('sratofastq/{srr}.fastq.gz', srr=ALL_SRR),
        dynamic_input = merge_fastq_input
    benchmark: 'benchmarks/merge_fastq/{sample}.txt'
    output: temp('fastq_merged/{sample}.fastq.gz')
    wrapper:
        get_wrapper('merge_fastq_wrapper')

rule quantify_cdna:
    input:
        R1='fastq_merged/{sample}.fastq.gz',
    output:
        'counts_cdna/{sample}/abundance.tsv',
    params:
        index=CDNA_IDX,
        outdir='counts_cdna/{sample}'
    threads: 16 ## hardcoded below
    shell:
        r'''
        kallisto quant \
        --index={params.index} \
        --threads={threads}\
        -l 200 \
        -s 30 \
        --output-dir={params.outdir} \
        --single \
        -b 100 <(zcat {input.R1})'''
rule quantify_cds:
    input:
        R1='fastq_merged/{sample}.fastq.gz',
    output:
        'counts_cds/{sample}/abundance.tsv',
    params:
        index=CDS_IDX,
        outdir='counts_cds/{sample}'
    threads: 16 ## hardcoded below
    shell:
        r'''
        kallisto quant \
        --index={params.index} \
        --threads={threads}\
        -l 200 \
        -s 30 \
        --output-dir={params.outdir} \
        --single \
        -b 100 <(zcat {input.R1})'''
