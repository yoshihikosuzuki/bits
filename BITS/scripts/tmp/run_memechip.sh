#!/bin/bash
#$ -o sge.log
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe smp 
#$ -l hostname=
#$ -N 

INPUT_PREFIX=$1   # hoge(.gff)
MOD_TYPE=$2   # [m6A|m4C|others|all]
THREAD_NUM=$3

## python ~/work2/pbmga/src/generate_fasta_for_memechip.py hoge.gff

meme-chip -oc ${INPUT_PREFIX}.memechip_${MOD_TYPE} -index-name ${INPUT_PREFIX}.memechip_${MOD_TYPE}.html -order 1 -meme-mod zoops -meme-minw 4 -meme-maxw 20 -meme-nmotifs 10 -meme-p ${THREAD_NUM} -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 ${INPUT_PREFIX}.memechip_${MOD_TYPE}.fasta
