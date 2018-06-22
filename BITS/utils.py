import sys
import subprocess
from logzero import logger


def run_command(command):
    """
    General-purpose command executer.
    """
    
    try:
        out = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as proc:
        logger.error(proc.output.decode('utf-8'))
        sys.exit(1)
    else:
        return out.decode('utf-8')
