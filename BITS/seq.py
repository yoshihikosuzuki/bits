import os
from logzero import logger
import matplotlib.pyplot as plt
import matplotlib.image as img

from .utils import run_command


def load_fasta(in_fname):
    from Bio.SeqIO import FastaIO
    with open(in_fname, 'r') as f:
        return dict(FastaIO.SimpleFastaParser(f))


RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def revcomp(seq):
    """
    Return the reverse complement of the given sequence.
    """

    return ''.join([RC_MAP[c] for c in seq[::-1]])


class DotPlot:
    def __init__(self, out_dir, gepard_command):
        if not os.path.isdir(out_dir):
            run_command(f"mkdir {out_dir}")

        self.out_dir = out_dir   # temporary directory
        self.gepard_command = gepard_command
        self.a_fname = os.path.join(out_dir, "a.fasta")
        self.b_fname = os.path.join(out_dir, "b.fasta")
        self.dotplot_fname = os.path.join(out_dir, "dotplot.png")

    def plot(self, a, b, a_name="a", b_name="b"):
        """
        Generate a dot plot of the given two sequences.
        """

        def write_fasta(seq, seq_name, fname):
            with open(fname, 'w') as f:
                f.write(f">{seq_name}/0/0_{len(seq)}\n{seq}\n")

        write_fasta(a, a_name, self.a_fname)
        write_fasta(b, b_name, self.b_fname)
    
        # Calculate dot plot
        command = f"unset DISPLAY; {self.gepard_command} -seq1 {self.a_fname} -seq2 {self.b_fname} -outfile {self.dotplot_fname}"
        run_command(command)
    
        # Show dot plot
        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.tick_params(labelleft=False, left=False)
        # this assignment and plt.show() are necessary to show only one figure
        plt.imshow(img.imread(self.dotplot_fname))
        plt.show()


def extract_adapters_from_bax(in_bax, out_adapters):
    """
    Extract all adapter sequences detected in PacBio reads
    """
    
    from pbcore.io import BasH5Reader

    with open(out_adapters, 'w') as out:
        with BasH5Reader(in_bax) as f:
            for r in f:
                for a in r.adapters:
                    out.write(a.basecalls() + '\n')


def consensus_adaptors(in_bax, out_adapters='adapters', out_consensus='adapter.consensus.fasta'):
    """
    Extract all adapter sequences from *.bax.h5, and then take consensus of them.
    """

    extract_adapters_from_bax(in_bax, out_adapters)
    #consed.consensus(out_adapters, out_consensus)   # TODO: update
