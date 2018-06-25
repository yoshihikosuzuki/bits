#!/bin/bash

HEADERS=preads.adaptor_unremoved
FASTA=preads4falcon.fasta

fatt extract --file ${HEADERS} ${FASTA} > ${HEADERS}.fasta
mummer -maxmatch -b -c -F ${HEADERS}.fasta ${HEADERS}.fasta > ${HEADERS}.mums
mummerplot -R ${HEADERS}.fasta -Q ${HEADERS}.fasta -postscript -p ${HEADERS} ${HEADERS}.mums
