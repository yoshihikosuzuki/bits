import sys
from .utils import slurm_nize

if __name__ == "__main__":
    with open(f"{sys.argv[1]}.slurm", 'w') as f:
        f.write(slurm_nize(f"bash {sys.argv[1]}\n",
                           **dict(arg.split('=') for arg in sys.argv[2:])))
