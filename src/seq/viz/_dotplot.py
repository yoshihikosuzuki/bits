from dataclasses import InitVar, dataclass, field
from os.path import join, splitext
from typing import List, Optional, Tuple, Type, Union

import plotly.graph_objects as go
from bits.util import run_command
from logzero import logger
from plotly_light import show_image

from .._io import FastaRecord, save_fasta


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
        assert gepard_root is not None or (
            gepard_jar is not None and gepard_mat is not None
        ), "`gepard_root` or `gepard_[jar|mat]` must be specified"
        _gepard_jar = (
            gepard_jar
            if gepard_jar is not None
            else f"{gepard_root}/dist/Gepard-1.40.jar"
        )
        _gepard_mat = (
            gepard_mat
            if gepard_mat is not None
            else f"{gepard_root}/resources/matrices/edna.mat"
        )
        self.gepard = " ".join(
            [
                f"java -cp {_gepard_jar}",
                "org.gepard.client.cmdline.CommandLine",
                f"-matrix {_gepard_mat}",
            ]
        )
        run_command(f"mkdir -p {self.tmp_dir}")

    def plot(
        self,
        a_seqs: Optional[Union[str, Type[FastaRecord], List[Type[FastaRecord]]]],
        b_seqs: Optional[Union[str, Type[FastaRecord], List[Type[FastaRecord]]]],
        a_range: Optional[Tuple[int, int]] = None,
        b_range: Optional[Tuple[int, int]] = None,
        word_size: int = 10,
        only_plot: bool = False,
        fig_size: Union[int, Tuple[Optional[int], Optional[int]]] = 1000,
        plot_size: int = 500,
        layout: go.Layout = None,
        out_fname: Optional[str] = None,
        static: bool = False,
    ):
        """Draw a dot plot between two sequences.

        positional arguments:
          @ [a|b]_seq : Sequence, FastaRecord, or fasta file name.

        optional arguments:
          @ word_size  : Word size for Gepard.
          @ only_plot  : Show only the plot and not texts and margins.
          @ fig_size   : Size of the png file of the dot plot (in pixel).
          @ plot_size  : Display size of the plot image (in pixel).
          @ out_fname  : Output file name of the dot plot.
        """

        def _prep(seqs: str, prolog: str):
            if (
                isinstance(seqs, str) and splitext(seqs)[1] == ".fasta"
            ):  # fasta file name
                return seqs

            out_fasta = join(self.tmp_dir, f"{prolog}.fasta")
            if not isinstance(seqs, list):
                seqs = [seqs]
            save_fasta(
                [
                    (
                        seq
                        if isinstance(seq, FastaRecord)
                        else FastaRecord(name=f"{prolog}/{i}/0_{len(seq)}", seq=seq)
                    )
                    for i, seq in enumerate(seqs)
                ],
                out_fasta,
            )
            return out_fasta

        a_fasta = _prep(a_seqs, "a")
        b_fasta = _prep(b_seqs, "b")

        if out_fname is None:
            out_fname = f"{self.tmp_dir}/dotplot.png"

        run_command(
            " ".join(
                [
                    "unset DISPLAY;",
                    f"{self.gepard}",
                    f"-seq1 {a_fasta}",
                    f"-seq2 {b_fasta}",
                    f"-maxwidth {fig_size}",
                    f"-maxheight {fig_size}",
                    f"-word {word_size}",
                    "-onlyplot" if only_plot else "",
                    f"-outfile {out_fname}",
                ]
            )
        )

        if not only_plot and (a_range is not None or b_range is not None):
            logger.info("`[a|b]_range` is ignored (used only when `only_plot`)")
            a_range, b_range = None, None

        show_image(
            out_fname,
            static=static,
            width=plot_size,
            height=plot_size,
            x_range=a_range,
            y_range=b_range,
            layout=layout,
        )
