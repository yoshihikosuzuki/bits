#!/bin/bash

N_ARGS=2
if [ $# -ne ${N_ARGS} ]; then
    echo '
Usage: ./add_jupyter_kernel.sh <venv_directory> <kernel_name>

Add a virtual environment as a kernel of Jupyter Notebook. To show the list of current Jupyter kernels, run `$ jupyter kernelspec list`.

Options
-------

venv_directory [str]
    Relative path to the directory of virtual environment.

kernel_name [str]
    The name of the kernel. Specify arbitrary, identifiable name.
'
    exit 1
fi

VENV_DIR=$1
KERNEL_NAME=$2

source ${VENV_DIR}/bin/activate
ipython kernel install --user --name=${KERNEL_NAME} --display-name=${KERNEL_NAME}
