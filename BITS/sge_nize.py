import sys
from .utils import sge_nize

if __name__ == "__main__":
    sge_nize(sys.argv[1], **dict(arg.split('=') for arg in sys.argv[2:]))
