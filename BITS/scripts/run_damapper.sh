#!/bin/bash
#$ -N damapper
#$ -o sge.log
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe smp 8
#$ -l hostname=ax02

fasta2DAM CONTIGS pilon.fasta
DBsplit -s500 CONTIGS
fasta2DB READS filtered_subreads.fasta
DBsplit -s500 READS
HPC.damapper -T8 -C -N CONTIGS READS | bash -v
LAcat CONTIGS.READS.*.las > CONTIGS.READS.las
LAsort -a CONTIGS.READS.las
LAdump -o -c CONTIGS READS CONTIGS.READS.S.las > ladump
DBdump -r -h CONTIGS > dbdump
