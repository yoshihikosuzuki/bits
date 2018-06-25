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

#mkdir metabat
#cd metabat
#ln -sf ../data/${REF} .

bowtie2-build -f ${REF} ${REF}.index

for DATA in ${SHORT_READ_DIR}/*.10m
do
    FNAME=$(echo ${DATA} | awk -F'/' '{print $NF}')
    bowtie2 --no-unal -p ${THREAD_NUM} -q -x ${REF}.index -U ${DATA} -S ${FNAME}.sam
    samtools view -@ ${THREAD_NUM} -bS ${FNAME}.sam > ${FNAME}.bam
    samtools sort -@ ${THREAD_NUM} ${FNAME}.bam -o ${FNAME}.sort.bam
done

jgi_summarize_bam_contig_depths --minMapQual 4 --outputDepth ${REF}.depth ./*.sort.bam
metabat --verysensitive --unbinned -i ${REF} -a ${REF}.depth -o ${REF%.fasta}.bin

python ${SCRIPT_DIR}extract_long_unbinned_contigs.py ${REF%.fasta}.bin.unbinned.fa