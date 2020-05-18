#!/bin/bash
# Template of a shell script

N_ARGS=
if [ $# -ne ${N_ARGS} ]; then
    echo '
Usage: ./script_name.sh <argument>

Description.

Options
-------

argument [type]
    Description.
'
    exit 1
fi

=$1


