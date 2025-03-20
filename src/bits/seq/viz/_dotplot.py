from dataclasses import InitVar, dataclass, field
from os.path import basename, join, splitext
from typing import List, Optional, Tuple, Type, Union

import plotly.graph_objects as go
import plotly_light as pl
from bits.util import run_command

from .._io import FastaRecord, load_fasta, save_fasta


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
        a_name: Optional[str] = None,
        b_name: Optional[str] = None,
        a_range: Optional[Tuple[int, int]] = None,
        b_range: Optional[Tuple[int, int]] = None,
        word_size: int = 10,
        fig_size: Union[int, Tuple[Optional[int], Optional[int]]] = 1000,
        plot_size: int = 500,
        show_grid: bool = False,
        grid_col: Optional[str] = None,
        opacity: Optional[float] = None,
        layout: go.Layout = None,
        return_fig: bool = False,
        out_fname: Optional[str] = None,
        original_plot: bool = False,
        static: bool = False,
    ):
        """Draw a dot plot between two sequences.

        positional arguments:
          @ [a|b]_seq : Sequence, FastaRecord, or fasta file name.

        optional arguments:
          @ word_size  : Word size for Gepard.
          @ fig_size   : Size of the png file of the dot plot (in pixel).
          @ plot_size  : Display size of the plot image (in pixel).
          @ out_fname  : Output file name of the dot plot.
          @ original_plot  : Show the original plot of Gepard.
          @ static     : Show the plot as a static image.
        """

        def _prep(seqs: str, prolog: str):
            # fasta file name
            if isinstance(seqs, str) and splitext(seqs)[1] in (
                ".fasta",
                ".fa",
            ):
                return seqs, basename(seqs)

            name = prolog
            if not isinstance(seqs, list):  # single sequence
                if isinstance(seqs, FastaRecord):
                    name = seqs.name
                seqs = [seqs]
            out_seqs = []
            for i, seq in enumerate(seqs):
                if isinstance(seq, FastaRecord):
                    out_seqs.append(seq)
                else:
                    out_seqs.append(
                        FastaRecord(name=f"{prolog}/{i}/0_{len(seq)}", seq=seq)
                    )
            out_fasta = join(self.tmp_dir, f"{prolog}.fasta")
            save_fasta(out_seqs, out_fasta, verbose=False)
            return out_fasta, name

        a_fasta, a_title = _prep(a_seqs, "a")
        b_fasta, b_title = _prep(b_seqs, "b")
        if a_name is not None:
            a_title = a_name
        if b_name is not None:
            b_title = b_name

        if original_plot:
            a_range = b_range = None
        else:
            if a_range is None:
                a_range = (
                    0,
                    sum([seq.length for seq in load_fasta(a_fasta, verbose=False)]),
                )
            if b_range is None:
                b_range = (
                    0,
                    sum([seq.length for seq in load_fasta(b_fasta, verbose=False)]),
                )
            # Flip `b_range` so that left top of the image is (0, 0)
            b_range = (b_range[1], b_range[0])

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
                    "-onlyplot" if not original_plot else "",
                    f"-outfile {out_fname}",
                ]
            )
        )

        fig = pl.show_image(
            out_fname,
            static=static,
            width=plot_size,
            height=plot_size,
            x_range=a_range,
            y_range=b_range,
            layer="above",
            opacity=opacity,
            layout=pl.merge_layout(
                pl.layout(
                    x_title=a_title,
                    y_title=b_title,
                    x_grid=show_grid,
                    y_grid=show_grid,
                    grid_col=grid_col,
                ),
                layout,
            ),
            return_fig=return_fig,
        )
        if return_fig:
            return fig

    def plot_self(
        self,
        seqs: Optional[Union[str, Type[FastaRecord], List[Type[FastaRecord]]]],
        name: Optional[str] = None,
        range: Optional[Tuple[int, int]] = None,
        word_size: int = 10,
        fig_size: Union[int, Tuple[Optional[int], Optional[int]]] = 1000,
        plot_size: int = 500,
        show_grid: bool = False,
        grid_col: Optional[str] = None,
        opacity: Optional[float] = None,
        layout: go.Layout = None,
        return_fig: bool = False,
        out_fname: Optional[str] = None,
        original_plot: bool = False,
        static: bool = False,
    ):
        """Draw a self-vs-self dot plot of a single sequence.

        positional arguments:
          @ seq : Sequence, FastaRecord, or fasta file name.

        optional arguments:
          @ word_size  : Word size for Gepard.
          @ fig_size   : Size of the png file of the dot plot (in pixel).
          @ plot_size  : Display size of the plot image (in pixel).
          @ out_fname  : Output file name of the dot plot.
          @ original_plot  : Show the original plot of Gepard.
          @ static     : Show the plot as a static image.
        """

        return self.plot(
            seqs,
            seqs,
            name,
            name,
            range,
            range,
            word_size,
            fig_size,
            plot_size,
            show_grid,
            grid_col,
            opacity,
            layout,
            return_fig,
            out_fname,
            original_plot,
            static,
        )
