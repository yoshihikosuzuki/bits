from os.path import join
from dataclasses import dataclass
import matplotlib.pyplot as plt
import matplotlib.image as img
from BITS.util.proc import run_command
from .io import save_fasta


@dataclass(repr=False, eq=False)
class DotPlot:
    """Draw dot plot between 2 strings or fasta files using Gepard."""
    gepard    : str           # path to the Gepard jar executable file
    fig_size  : int = 750
    plot_size : int = 11
    word_size : int = 10
    out_dir   : str = "tmp"   # directory to temporarily put fasta and png files

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")
        self.out_fname = join(self.out_dir, "dotplot.png")

    def _plot(self, a_fname, b_fname):
        run_command(' '.join([f"unset DISPLAY;",
                              f"{self.gepard} -seq1 {a_fname} -seq2 {b_fname}",
                              f"-maxwidth {self.fig_size} -maxheight {self.fig_size}",
                              f"-word {self.word_size}",
                              f"-outfile {self.out_fname}"]))
        fig, ax = plt.subplots(figsize=(self.plot_size, self.plot_size))
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
