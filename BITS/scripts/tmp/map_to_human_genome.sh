#!/bin/bash

HUMANREF=$1
SUBREADS=$2
THREAD_NUM=$3   # MUST be power of 2
ENV_ACTICATE=$4

source ${ENV_ACTICATE}
fasta2DAM ref ${HUMANREF}
fasta2DB subreads ${SUBREADS}
DBsplit -x500 -s500 ref
DBsplit -x500 -s500 subreads
HPC.damapper -T${THREAD_NUM} ref subreads | bash -v
LAcat subreads.#.ref.las > subreads.ref.las
LA4Falcon -m subreads ref subreads.ref > ovl.all
