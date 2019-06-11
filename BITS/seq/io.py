from Bio.SeqIO import FastaIO
from .utils import split_seq


def load_fasta(in_fname):
    """Load a fasta file as a dict of {header: seq}."""
    with open(in_fname, 'r') as f:
        return dict(FastaIO.SimpleFastaParser(f))


def save_fasta(reads, out_fname, sort=True, width=-1):
    """
    <reads> must be a dict of {header: seq}.
    If <sort> is True, the headers will be sorted.
    Newlines are inserted at every <width> bp in each sequence (-1 means no newlines).
    """
    with open(out_fname, 'w') as f:
        for header, seq in sorted(reads.items()) if sort else reads.items():
            if width > 0:
                seq = '\n'.join(split_seq(seq, width))
            f.write(f">{header}\n{seq}\n")
