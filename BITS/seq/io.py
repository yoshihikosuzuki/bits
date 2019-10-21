import numpy as np
from Bio import SeqIO
import dnaio
from .utils import ascii_to_phred, split_seq


def load_fasta(in_fname):
    """Load a fasta file as a dictionary of {header: seq}."""
    with open(in_fname, 'r') as f:
        return dict(SeqIO.FastaIO.SimpleFastaParser(f))


def load_fastq(in_fname, only_qual=False):
    """Load a fastq file as a dictionary of `{header: (seq, qual)}`, where `qual` is positional QVs.
    If <only_qual> is True, it instead returns {header: qual}."""
    seqs = {}
    with dnaio.open(in_fname) as f:
        for record in f:
            quals = np.array([ascii_to_phred(a) for a in record.qualities])
            seqs[record.name] = quals if only_qual else (record.sequence, quals)
    return seqs


def save_fasta(reads, out_fname, sort=True, width=-1):
    """<reads> must be a dict of {header: seq}.
    If <sort> is True, the headers will be sorted.
    Newlines are inserted at every <width> bp in each sequence (-1 means no newlines)."""
    with open(out_fname, 'w') as f:
        for header, seq in sorted(reads.items()) if sort else reads.items():
            if width > 0:
                seq = '\n'.join(split_seq(seq, width))
            f.write(f">{header}\n{seq}\n")
