#!/bin/bash
#$ -N JOB_NAME
#$ -o OUT_LOG
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe smp N_CORE

REF_FASTA=$1
READS_FASTA=$2

fasta2DAM REF ${REF_FASTA}
DBsplit -s500 REF
fasta2DB READS ${READS_FASTA}
DBsplit -s500 READS
HPC.damapper -T8 -C -N REF READS | bash -v
LAcat REF.READS.*.las > REF.READS.las
LAsort -a REF.READS.las
LAdump -o -c REF READS REF.READS.S.las > REF.READS.ladump
DBdump -r -h REF > REF.dbdump
