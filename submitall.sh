#!/bin/bash
source activate clipseq
echo $PWD

snakemake --js $PWD/jobscript.sh\
    --printshellcmds\
    --cluster-config $PWD/cluster.yaml\
    --jobname 'make.{jobid}.{rulename}'\
    --keep-going\
    --stats $PWD/snakemake.stats\
    --timestamp\
    --rerun-incomplete\
    -j 100\
    --cluster 'qsub -q cmb -l walltime={cluster.time} -l mem={cluster.mem} -l vmem={cluster.mem} -l pmem={cluster.mem} -l nodes=1:ppn={cluster.cores} -o {cluster.logdir} -e {cluster.logdir}'
