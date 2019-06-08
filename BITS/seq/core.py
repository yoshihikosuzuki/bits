from os.path import join
from dataclasses import dataclass
from logzero import logger
import matplotlib.pyplot as plt
import matplotlib.image as img
from BITS.util.proc import run_command

RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def single_to_multi(seq, width=100):
    """
    Cut a single sequence at every <width> bp and return as a list.
    """

    return [seq[i:i + width] for i in range(0, len(seq), width)]


def revcomp(seq):
    return ''.join([RC_MAP[c] for c in seq[::-1]])


def compress_homopolymer(seq):
    ret = ""
    prev = ""
    for s in seq:
        if s != prev:
            ret += s
        prev = s
    return ret


@dataclass(repr=False)
class Cigar:
    string: str

    def iter(self):
        """
        Generator of a series of tuples (num. of bases, operator).
        """

        length = ""
        for c in self.string:
            if c == '=' or c == 'D' or c == 'I' or c == 'X' or c == 'N':
                yield (int(length), c)
                length = ""
            else:
                length += c

    @property
    def alignment_len(self):
        """
        Alignment length including gaps and masked regions.
        """

        return sum([l for l, op in self.iter()])

    def reversed(self):
        """
        Just reverse CIGAR without swapping in/del. Used for reverse complement.
        """

        return ''.join(reversed([f"{l}{op}" for l, op in self.iter()]))

    def flatten(self):
        """
        Convert CIGAR to a sequence of operations. Used for masking some specific intervals.
        """

        return FlattenCigar(''.join([op for l, op in self.iter() for i in range(l)]))

    def mask_intvl(self, intvl, ignore_op='D'):
        # TODO: how to do about I/D/X around intervals' boundaries

        cigar_f = self.flatten()
    
        starts = [i[0] for i in intvl]
        ends = [i[1] for i in intvl]
        index = 0
        pos = 0
        for i, c in enumerate(cigar_f):
            if index >= len(starts):
                break
            if i != 0 and c != ignore_op:
                pos += 1
            if pos > ends[index]:
                index += 1
            if index < len(starts) and starts[index] <= pos and pos <= ends[index]:   # NOTE: end-inclusive
                cigar_f[i] = 'N'

        return cigar_f.unflatten()


@dataclass(repr=False)
class FlattenCigar:
    string: str

    def unflatten(self):
        """
        Convert a series of oprations to a normal CIGAR string.
        """

        cigar = ""
        count = 0
        prev_c = self.string[0]
        for c in self.string:
            if c == prev_c:
                count += 1
            else:
                cigar += f"{count}{prev_c}"
                count = 1
                prev_c = c
        cigar += f"{count}{prev_c}"
        return Cigar(cigar)


@dataclass(repr=False, eq=False)
class DotPlot:
    out_dir: str
    gepard_command: str

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")
        self.dotplot_fname = join(self.out_dir, "dotplot.png")

    def _plot(self, a_fname, b_fname):
        run_command(' '.join([f"unset DISPLAY;",
                              f"{self.gepard_command}",
                              f"-seq1 {a_fname}",
                              f"-seq2 {b_fname}",
                              f"-outfile {self.dotplot_fname}"]))

        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.tick_params(labelleft=False, left=False)
        plt.imshow(img.imread(self.dotplot_fname))
        plt.show()

    def plot_fasta(self, a_fname, b_fname):
        self._plot(a_fname, b_fname)

    def plot(self, a_seq, b_seq, a_name="a", b_name="b"):
        a_fname, b_fname = join(self.out_dir, "a.fasta"), join(self.out_dir, "b.fasta")
        save_fasta({f"{a_name}/0/0_{len(a_seq)}": a_seq}, a_fname)
        save_fasta({f"{b_name}/0/0_{len(b_seq)}": b_seq}, b_fname)
        self._plot(a_fname, b_fname)
