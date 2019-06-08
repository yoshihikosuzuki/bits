from .core import single_to_multi


def load_fasta(in_fname):
    """
    Load a fasta file as a dict of {header: seq}.
    """

    from Bio.SeqIO import FastaIO
    with open(in_fname, 'r') as f:
        return dict(FastaIO.SimpleFastaParser(f))


def save_fasta(reads, out_fname, sort=True, out_type="single", width=100):
    """
    NOTE: <reads> must be a dictionary.
    If <sort> is True, the headers in <reads> will be sorted.
    <out_type> defines the existence of newlines within the sequences (by every <width> bp).
    """

    assert out_type in set(["single", "multi"]), "<out_type> must be 'single' or 'multi'."
    with open(out_fname, 'w') as f:
        for header, seq in sorted(reads.items()) if sort else reads.items():
            if out_type == "multi":
                seq = '\n'.join(single_to_multi(seq, width))
            f.write(f">{header}\n{seq}\n")
