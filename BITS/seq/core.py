from os.path import join
from dataclasses import dataclass
import matplotlib.pyplot as plt
import matplotlib.image as img
from BITS.util.proc import run_command
from .io import save_fasta

RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))
CIGAR_CHAR = set(['=', 'D', 'I', 'X', 'N'])


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

    def __str__(self):
        return self.string

    def __iter__(self):
        self._objs = []   # list of the tuples (count, cigar) in <self.string>
        count = ""
        for c in self.string:
            if c in CIGAR_CHAR:
                self._objs.append((int(count), c))
                count = ""
            else:
                count += c
        self._i = 0   # index of <self._objs>
        return self

    def __next__(self):
        if self._i == len(self._objs):
            raise StopIteration()
        ret = self._objs[self._i]
        self._i += 1
        return ret

    @property
    def alignment_len(self):
        """Alignment length including gaps and masked regions."""
        return sum([l for l, c in self])

    def reverse(self):
        """Just reverse CIGAR without swapping in/del. Used for reverse complement."""
        return ''.join(reversed([f"{l}{c}" for l, c in self]))

    def swap_indel(self):
        """Swap I and D. This inverts the role of query and target."""
        self.string = self.string.replace("I", "?").replace("D", "I").replace("?", "D")

    def flatten(self):
        """Convert to a sequence of the operations."""
        return FlattenCigar(''.join([c for l, c in self for i in range(l)]))

    def mask_intvl(self, intvl, ignore_op='D'):   # TODO: refactor
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
    """Class for representing CIGAR as a sequence of operations, which is easier to handle."""
    string: str

    def __str__(self):
        return self.string

    def __iter__(self):
        self._i = 0
        return self

    def __next__(self):
        if self._i == len(self.string):
            raise StopIteration()
        ret = self.string[self._i]
        self._i += 1
        return ret

    def unflatten(self):
        """Convert to the normal CIGAR string."""
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
    """Draw dot plot between 2 strings/fasta files using Gepard."""
    gepard: str
    out_dir: str = "tmp"

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")
        self.out_fname = join(self.out_dir, "dotplot.png")

    def _plot(self, a_fname, b_fname):
        run_command(' '.join([f"unset DISPLAY;",
                              f"{self.gepard} -seq1 {a_fname} -seq2 {b_fname}",
                              f"-outfile {self.out_fname}"]))
        fig, ax = plt.subplots(figsize=(11, 11))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.tick_params(labelleft=False, left=False)
        plt.imshow(img.imread(self.out_fname))
        plt.show()

    def plot_fasta(self, a_fname, b_fname):
        """Show dot plot between 2 fasta files."""
        self._plot(a_fname, b_fname)

    def plot(self, a_seq, b_seq, a_name="a", b_name="b"):
        """Show dot plot between 2 strings."""
        a_fname, b_fname = join(self.out_dir, "a.fasta"), join(self.out_dir, "b.fasta")
        save_fasta({f"{a_name}/0/0_{len(a_seq)}": a_seq}, a_fname)
        save_fasta({f"{b_name}/0/0_{len(b_seq)}": b_seq}, b_fname)
        self._plot(a_fname, b_fname)
