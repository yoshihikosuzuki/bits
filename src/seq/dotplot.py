from os.path import join
from dataclasses import dataclass, field, InitVar
from typing import Type, Union, Optional
from BITS.plot.matplotlib import show_image
from BITS.util.proc import run_command
from .io import FastaRecord, save_fasta


@dataclass(repr=False, eq=False)
class DotPlot:
    """Draw dot plot between two sequences or fasta files with Gepard.

    usage in Jupyter:
      > gepard_jar = "/path/to/gepard/dist/Gepard-1.40.jar"
      > gepard_mat = "/path/to/gepard/resources/matrices/edna.mat"
      > dp = DotPlot(gepard_jar, gepard_mat)
      > dp.plot(seq1, seq2)

    positional variables:
      @ gepard_jar : Path to a jar executable of Gepard.
      @ gepard_mat : Path to a score matrix of Gepard.

    optional variables:
      @ tmp_dir : Directory for generating fasta files and plots.
    """
    gepard_jar: InitVar[str]
    gepard_mat: InitVar[str]
    gepard: str = field(init=False)
    tmp_dir: str = "tmp"

    def __post_init__(self, gepard_jar, gepard_mat):
        run_command(f"mkdir -p {self.tmp_dir}")
        self.gepard = ' '.join([f"java -cp {gepard_jar}",
                                "org.gepard.client.cmdline.CommandLine",
                                f"-matrix {gepard_mat}"])

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
        show_image(out_fname, plot_size, plot_size)

    def plot(self,
             a_seq: Union[str, Type[FastaRecord]],
             b_seq: Union[str, Type[FastaRecord]],
             from_fasta: bool = False,
             a_name: Optional[str] = None,
             b_name: Optional[str] = None,
             out_fname: Optional[str] = None,
             word_size: int = 10,
             fig_size: int = 750,
             plot_size: int = 11):
        """Draw a dot plot between two sequences.

        positional arguments:
          @ [a|b]_seq : Sequence, FastaRecord, or fasta file name.

        optional arguments:
          @ from_fasta : `[a|b]_seq` are fasta file names.
          @ [a|b]_name : Display names of the sequences in the plot.
          @ out_fname  : Output file name of the dot plot.
          @ word_size  : Word size for Gepard.
          @ fig_size   : Size of the png file of the dot plot.
          @ plot_size  : Display size of the plot in Jupyter.
        """
        def _prep(seq: str, name: str, prolog: str):
            save_fasta([FastaRecord(name=(name if name is not None
                                          else seq.name if isinstance(seq, FastaRecord)
                                          else f"{prolog}/0/0_{len(seq)}"),
                                    seq=(seq.seq if isinstance(seq, FastaRecord)
                                         else seq))],
                       join(self.tmp_dir, f"{prolog}.fasta"))

        if not from_fasta:
            _prep(a_seq, a_name, 'a')
            _prep(b_seq, b_name, 'b')
        self._plot(a_seq if from_fasta else join(self.tmp_dir, "a.fasta"),
                   b_seq if from_fasta else join(self.tmp_dir, "b.fasta"),
                   (out_fname if out_fname is not None
                    else join(self.tmp_dir, "dotplot.png")),
                   word_size, fig_size, plot_size)
