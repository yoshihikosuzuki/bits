#!/bin/bash

source ~/.bash_profile

NAME=000022F.gff

cp data/motif_analysis/${NAME} ./
python ~/work2/pbmga/src/generate_fasta_for_memechip.py ${NAME}
~/work2/pbmga/src/shell_scripts/run_memechip.sh ${NAME%.gff} m4C 12
