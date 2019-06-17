import pandas as pd
from Bio.SeqIO import FastaIO
from BITS.util.proc import run_command
from .utils import split_seq


def load_fasta(in_fname):
    """Load a fasta file as a pd.DataFrame of [header: seq]."""
    with open(in_fname, 'r') as f:
        ret = pd.DataFrame.from_dict(dict(FastaIO.SimpleFastaParser(f)),
                                     orient="index", columns=["seq"])
        ret.index.names = ["header"]
        return ret


def save_fasta(reads, out_fname, sort=True, width=-1):
    """<reads> must be a dict of {header: seq}.
    If <sort> is True, the headers will be sorted.
    Newlines are inserted at every <width> bp in each sequence (-1 means no newlines).
    """
    with open(out_fname, 'w') as f:
        for header, seq in sorted(reads.items()) if sort else reads.items():
            if width > 0:
                seq = '\n'.join(split_seq(seq, width))
            f.write(f">{header}\n{seq}\n")


def load_dazz_db(in_fname):
    """Load DAZZ_DB file as a pd.DataFrame of [dbid: header seq]. Headers must not contain tab."""
    command = (f"DBshow {in_fname} | "
               f"awk 'BEGIN {{first = 1}} "
               f"{{if (substr($0, 1, 1) == \">\") "
               f"{{if (first == 1) {{first = 0}} else {{printf(\"%s\\t%s\\n\", header, seq)}}; "
               f"header = substr($0, 2); seq = "";}} "
               f"else {{seq = seq $0}}}} "
               f"END {{printf(\"%s\\t%s\\n\", header, seq)}}'")
    ret = pd.DataFrame.from_dict({dbid: tuple(*line.split('\t'))
                                  for dbid, line
                                  in enumerate(run_command(command).strip().split('\n'), start=1)},
                                 orient="index", columns=["header", "seq"])
    ret.index.names = ["dbid"]
    return ret
