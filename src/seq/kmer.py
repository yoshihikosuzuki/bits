from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from logzero import logger
from .util import revcomp_seq
from BITS.plot.plotly import make_scatter, make_layout, show_plot
from BITS.util.proc import run_command


def canonical_kmer(kmer):
    """Convert into a canonical k-mer."""
    return min(kmer, revcomp_seq(kmer))


def kmer_spectrum(seq, k, by="forward"):
    """Return a set of forward/canonical k-mers contained in `seq`.

    positional arguments:
      @ seq : Input sequence.
      @ k   : Length of k-mers.

    optional arguments:
      @ by  : Must be one of {"forward", "canonical"}.
    """
    return set([seq[i:i + k] if by == "forward"
                else canonical_kmer(seq[i:i + k])
                for i in range(len(seq) - k + 1)])


def create_meryl_db(fasta_fname, k, meryl_out):
    """Make a meryl k-mer database from the sequences in `fasta_fname`."""
    # TODO: `count` command of meryl counts canonical k-mers.
    #       For forward k-mers, use `count-forward`.
    run_command(f"meryl count k={k} {fasta_fname} output {meryl_out}")


def count_kmers(fasta_fname, k, greater_than=0, lower=True, out_dir=".", out_fname=None):
    """Count `k`-mer contained in `fasta_fname` using meryl. Meryl must be installed.
    NOTE: large `k` value might cause segmentation fault of meryl.
    NOTE: Meryl uses approximate counting. That is, some k-mers will not be counted.
    """
    assert isinstance(k, int), "k must be integer"
    if k % 2 == 0:
        logger.warn("Even value of `k` results in non-canonical k-mers")

    meryl_out = f"{out_dir}/{fasta_fname}.meryl" if out_fname is None else f"{out_dir}/{out_fname}"
    create_meryl_db(fasta_fname, k, meryl_out)

    kmer_counts = Counter()
    ret = run_command(f"meryl print greater-than {greater_than} {meryl_out}")
    for line in ret.strip().split('\n'):
        kmer, count = line.split('\t')
        kmer_counts[kmer.lower() if lower else kmer] = int(count)
    return kmer_counts


def hist_kmer_counts(kmer_counts, min_count=1, max_count=100, log_scale=True):
    """X-axis is count of a k-mer, y-axis is frequency of k-mers having the count."""
    plt.hist(list(kmer_counts.values()), bins=max_count - min_count + 1,
             range=(min_count, max_count), log=log_scale)
    plt.grid(b=True, which="both")
    plt.show()


def seq_to_kmer_count_spectrum(seq, kmer_counts):
    """Convert a DNA sequence into a spectrum of positional k-mer counts, which represents
    repetitiveness and the degree of error of the k-mers, i.e., the positions."""
    k = len(next(iter(kmer_counts.keys())))
    return np.array([kmer_counts[canonical_kmer(seq[i:i + k])]
                     for i in range(len(seq) - k + 1)])


def plot_kmer_count_spectrum(seq, kmer_counts, y_max=None):
    """Plot the k-mer spectrum of `seq` according to `kmer_counts`."""
    spectrum = seq_to_kmer_count_spectrum(seq, kmer_counts)
    show_plot([make_scatter(np.arange(len(spectrum)), spectrum,
                            marker_size=3, show_legend=False)],
              make_layout(y_range=None if y_max is None else (0, y_max)))


def plot_kmer_count_spectrums(seqs, kmer_counts, y_max=None):
    """Plot the k-mer spectrum of `seqs` according to `kmer_counts`."""
    spectrums = [seq_to_kmer_count_spectrum(seq, kmer_counts) for seq in seqs]
    show_plot([make_scatter(np.arange(len(spectrum)), spectrum,
                            marker_size=3, name=f"Seq {i}")
               for i, spectrum in enumerate(spectrums)],
              make_layout(x_range=(0, max([len(spectrum) for spectrum in spectrums])),
                          y_range=None if y_max is None else (0, y_max)))
