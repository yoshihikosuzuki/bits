#!/bin/bash
#$ -o sge.log
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe smp 
#$ -l hostname=
#$ -N 

THREAD_NUM=
#REF=pilon.fasta


#SHORT_READ_DIR=/work2/yoshihiko_s/metaFALCON/paper/japanese_10m
#SCRIPT_DIR=/work2/yoshihiko_s/metaFALCON/paper/assembly/unitig/


source ~/.bash_profile
#source ~/work2/env1/bin/activate

#rm -r ${REF}.checkm checkm_result
checkm lineage_wf -t ${THREAD_NUM} -x fa -f ${REF}.checkm . ./checkm_result