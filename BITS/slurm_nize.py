import sys
from .utils import slurm_nize

if __name__ == "__main__":
    slurm_nize(sys.argv[1], **dict(arg.split('=') for arg in sys.argv[2:]))
