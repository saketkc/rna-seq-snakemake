# rna-seq-snakamake
Snakemake based pipeline for RNA-Seq analysis

Inputs to all snakefiles are specified using a `config.py`, and generally is as simple as specifiying the directory
containing all `*.fq/*.fq.gz/*.fastq.gz/*.sra` files.

1. [Snakefile_fastq_multilane.snake](Snakefile_fastq_multilane.snake) - For multilane paired end fastq files

2. [Snakefile_sra_pe](Snakefile_sra_pe) - For paired end `*.sra` files

3. [Snakefile_sra_se](Snakefile_sra_se) - For single end `*.sra` files

4. [Snakefile_fastq](Snakefile_fastq) - For paired end fastq files

## Dependecies

- snakemake
- STAR
- bamtools
- DESeq2
- qualimap
- multiqc
- rseqc
