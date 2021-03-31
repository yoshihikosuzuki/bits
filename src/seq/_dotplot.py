from os.path import join, splitext
from dataclasses import dataclass, field, InitVar
from typing import Type, Union, Optional
from ..util import run_command, show_image
from ._io import FastaRecord, save_fasta


@dataclass(repr=False, eq=False)
class DotPlot:
    """Draw dot plot between two sequences or fasta files with Gepard.

    Usage example in Jupyter:
      > gepard_root = "/path/to/gepard"
      > dp = DotPlot(gepard_root)
      > dp.plot(seq1, seq2)

    positional arguments:
      (One of `gepard_root` or `gepard_[jar|mat]` must be specified.)
      @ gepard_root : Path to the root directory of Gepard.
                      If specified, `f"{gepard_root}/dist/Gepard-1.40.jar"` and
                      `f"{gepard_root}/resources/matrices/edna.mat"` are used as
                      `gepard_jar` and `gepard_mat`.
      @ gepard_jar  : Path to a jar executable of Gepard.
      @ gepard_mat  : Path to a score matrix of Gepard.

    optional arguments:
      @ tmp_dir : Directory for generating fasta files and plots.
    """
    gepard_root: InitVar[Optional[str]] = None
    gepard_jar: InitVar[Optional[str]] = None
    gepard_mat: InitVar[Optional[str]] = None
    gepard: str = field(init=False)
    tmp_dir: str = "tmp"

    def __post_init__(self, gepard_root, gepard_jar, gepard_mat):
        assert gepard_root is not None or (gepard_jar is not None and gepard_mat is not None), \
            "`gepard_root` or `gepard_[jar|mat]` must be specified"
        _gepard_jar = (gepard_jar if gepard_jar is not None
                       else f"{gepard_root}/dist/Gepard-1.40.jar")
        _gepard_mat = (gepard_mat if gepard_mat is not None
                       else f"{gepard_root}/resources/matrices/edna.mat")
        self.gepard = ' '.join([f"java -cp {_gepard_jar}",
                                "org.gepard.client.cmdline.CommandLine",
                                f"-matrix {_gepard_mat}"])
        run_command(f"mkdir -p {self.tmp_dir}")

    def _plot(self,
              a_fname: str,
              b_fname: str,
              out_fname: str,
              word_size: int,
              fig_size: int,
              plot_size: int):
        run_command(' '.join(["unset DISPLAY;",
                              f"{self.gepard}",
                              f"-seq1 {a_fname}",
                              f"-seq2 {b_fname}",
                              f"-maxwidth {fig_size}",
                              f"-maxheight {fig_size}",
                              f"-word {word_size}",
                              f"-outfile {out_fname}"]))
        show_image(out_fname, plot_size)

    def plot(self,
             a_seq: Union[str, Type[FastaRecord]],
             b_seq: Union[str, Type[FastaRecord]],
             out_fname: Optional[str] = None,
             word_size: int = 10,
             fig_size: int = 750,
             plot_size: int = 11):
        """Draw a dot plot between two sequences.

        positional arguments:
          @ [a|b]_seq : Sequence, FastaRecord, or fasta file name.

        optional arguments:
          @ out_fname  : Output file name of the dot plot.
          @ word_size  : Word size for Gepard.
          @ fig_size   : Size of the png file of the dot plot (in pixel).
          @ plot_size  : Display size of the plot image (in inch).
        """
        def _prep(seq: str,
                  prolog: str):
            if isinstance(seq, str) and splitext(seq)[1] == ".fasta":
                return seq
            out_fname = join(self.tmp_dir, f"{prolog}.fasta")
            save_fasta(seq if isinstance(seq, FastaRecord)
                       else FastaRecord(name=f"{prolog}/0/0_{len(seq)}",
                                        seq=seq),
                       out_fname)
            return out_fname

        self._plot(_prep(a_seq, 'a'),
                   _prep(b_seq, 'b'),
                   (out_fname if out_fname is not None
                    else join(self.tmp_dir, "dotplot.png")),
                   word_size,
                   fig_size,
                   plot_size)
