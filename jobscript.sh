#!/bin/bash
# properties = {properties}

cd $PBS_O_WORKDIR
export PATH="/home/cmb-panasas2/skchoudh/software_frozen/anaconda2/envs/clipseq/bin:/home/cmb-panasas2/skchoudh/software_frozen/meme_4.10.2_bin/bin:/home/cmb-panasas2/skchoudh/software_frozen/anaconda2/bin:/home/cmb-panasas2/skchoudh/software_frozen/phast-1.4/bin:/home/cmb-panasas2/skchoudh/software_frozen/meme_4.10.2_bin/bin:/home/cmb-panasas2/skchoudh/software_frozen/anaconda2/bin:/home/cmb-panasas2/skchoudh/software_frozen/phast-1.4/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin"
source activate clipseq

echo $JOB_ID
echo "=============================="

{exec_job}

# Report resource consumption because it's not reported by default
echo "------------------------------"
qstat -f $JOB_ID | grep 'used'

# if the job succeeds, snakemake 
# touches jobfinished, thus if it exists cat succeeds. if cat fails, the error code indicates job failure
# an error code of 100 is needed since UGER only prevents execution of dependent jobs if the preceding
# job exits with error code 100

cat $1 &>/dev/null && exit 0 || exit 100
