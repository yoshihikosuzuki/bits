#!/bin/bash

SCRIPT_FILE=$1

# When using directly from shell
unalias qsub 2> /dev/null
unalias qsuball 2> /dev/null

source ~/.bash_profile
unset module   # suppress warnings appearing in the SGE log file (maybe due to above unaliasing?)
qsub ${SCRIPT_FILE}
