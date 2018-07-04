#!/bin/bash
#$ -N JOB_NAME
#$ -o OUT_LOG
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -pe smp N_CORE

sleep 24h
