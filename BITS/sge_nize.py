import sys
from .utils import sge_nize

if __name__ == "__main__":
    with open(f"{sys.argv[1]}.sge", 'w') as f:
        f.write(sge_nize(f"bash {sys.argv[1]}\n",
                         **dict(arg.split('=') for arg in sys.argv[2:])))
