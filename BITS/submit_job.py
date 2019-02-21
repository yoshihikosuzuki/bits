import sys
from .utils import submit_job

if __name__ == "__main__":
    submit_job(f"bash {sys.argv[1]}\n",
               f"{sys.argv[1]}.{sys.argv[2]}",
               sys.argv[2],
               sys.argv[3],
               **dict(arg.split('=') for arg in sys.argv[4:]))
