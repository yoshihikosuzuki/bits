import os.path.join
from dataclasses import dataclass
from logzero import logger
import matplotlib.pyplot as plt
import matplotlib.image as img
from .utils import run_command


def load_fasta(in_fname):
    """
    Load a fasta file as a dictionary of {header: seq}.
    """

    from Bio.SeqIO import FastaIO
    with open(in_fname, 'r') as f:
        return dict(FastaIO.SimpleFastaParser(f))


def single_to_multi(s, width=100):
    """
    Cut a single sequence at every <width> bp and return as a list.
    """

    return [s[i:i + width] for i in range(0, len(s), width)]


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


RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))

def revcomp(seq):
    return ''.join([RC_MAP[c] for c in seq[::-1]])


def homopolymer_compression(seq):
    ret = ""
    prev = ""
    for s in seq:
        if s != prev:
            ret += s
        prev = s
    return ret


@dataclass(repr=False, eq=False)
class DotPlot:
    out_dir: str
    gepard_command: str

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")
        self.a_fname = os.path.join(self.out_dir, "a.fasta")
        self.b_fname = os.path.join(self.out_dir, "b.fasta")
        self.dotplot_fname = os.path.join(self.out_dir, "dotplot.png")

    def plot(self, a, b, a_name="a", b_name="b"):
        """
        Generate a dot plot of the given two sequences.
        """

        save_fasta({f"{a_name}/0/0_{len(a)}": a}, self.a_fname)
        save_fasta({f"{b_name}/0/0_{len(b)}": b}, self.b_fname)

        # Calculate dot plot
        run_command(' '.join([f"unset DISPLAY;",
                              f"{self.gepard_command}",
                              f"-seq1 {self.a_fname}",
                              f"-seq2 {self.b_fname}",
                              f"-outfile {self.dotplot_fname}"]))
    
        # Show dot plot
        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.tick_params(labelleft=False, left=False)
        plt.imshow(img.imread(self.dotplot_fname))
        plt.show()
