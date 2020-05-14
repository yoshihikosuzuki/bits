#!/bin/bash

N_ARGS=2
if [ $# -ne ${N_ARGS} ]; then
    echo '
Usage: ./delete_jupyter_kernel.sh <kernel_name>

Delete a kernel of Jupyter Notebook. To show the list of current Jupyter kernels, run `$ jupyter kernelspec list`.

Options
-------

kernel_name [str]
    The name of a kernel to be deleted.
'
    exit 1
fi

KERNEL_NAME=$1

jupyter kernelspec uninstall ${KERNEL_NAME}
