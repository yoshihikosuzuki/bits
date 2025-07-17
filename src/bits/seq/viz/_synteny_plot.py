from typing import Callable, Dict, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly_light as pl
from pygenomeviz import GenomeViz

from .._type import BedRecord, SeqRecord


def load_syri(syri_file, min_len_indel=1000, min_len_sr=10000):
    """Load syri's output whose types are in `valid_types` and lengths are greater
    than `min_len_indel` (for allelic insertions and deletions) or `min_len_sr` (for
    the other distant structural rearrangements).

    Parameters
    ----------
    min_len_indel
        Minimum required length for INS and DEL (same allele).

    min_len_sr
        Minimum required length for the other structural rearrangements
        (different allele).
    """
    col_dtypes = {
        "RefChr": str,
        "RefStart": "Int64",
        "RefEnd": "Int64",
        "QryChr": str,
        "QryStart": "Int64",
        "QryEnd": "Int64",
        "Type": str,
    }
    df = pd.read_csv(
        syri_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 5, 6, 7, 10],
        names=list(col_dtypes.keys()),
        dtype=col_dtypes,
        na_values=["-"],
    )

    # filter by type
    # NOTE: "HDR" (highly divergent regions), "TDM" (tandem repeats) could also be
    #       used as annotations on the genomes
    valid_types = {"SYN", "INS", "DEL", "TRANS", "DUP", "INV", "INVTR", "INVDP"}
    df = df[df["Type"].isin(valid_types)]

    # filter by length
    df = df[
        (
            df["Type"].isin(["INS", "DEL"])
            & (
                np.maximum(df["RefEnd"] - df["RefStart"], df["QryEnd"] - df["QryStart"])
                >= min_len_indel
            )
        )
        | (
            ~df["Type"].isin(["INS", "DEL"])
            & (
                np.maximum(df["RefEnd"] - df["RefStart"], df["QryEnd"] - df["QryStart"])
                >= min_len_sr
            )
        )
    ]

    df = df.reset_index()
    return df


class SyntenyPlotGV(GenomeViz):
    """Wrapper of pyGenomeViz for genome plot.

    Example usage:
    ```
    > seq_names = ["A", "B", "C"]
    > seqs = [bs.load_fasta(f"name.fasta")[0] for name in seq_names]
    > anno_by_name = {name: bs.load_bed("anno.bed") for name in seq_names})
    > aln_by_pair = {("A", "B"): load_syri("A_vs_B.syri.out"), ("B", "C"): load_syri("B_vs_C.syri.out")}
    > cov = {name: bs.load_bed(f"{name}.cov.bedgraph", attrs=[("cov", float)]) for name in seq_names}
    >
    > gv = bs.SyntenyPlot(seqs, seq_names, show_scale_xticks=True)
    > gv.add_feature(anno_by_name)
    > gv.add_aln(aln_by_pair)
    > gv.add_subtrack(cov_by_name)
    > gv.show()
    ```
    """

    color_map = {
        "SYN": "gray",
        "TRANS": pl.colors["green"],
        "DUP": pl.colors["blue"],
        "INS": pl.colors["white"],
        "DEL": pl.colors["white"],
        "INV": pl.colors["yellow"],
        "INVTR": pl.colors["darkyellow"],
        "INVDP": pl.colors["darkyellow"],
    }

    def __init__(
        self,
        seqs: Sequence[SeqRecord],
        seq_names: Optional[Sequence[str]] = None,
        fig_width: float = 10,
        fig_track_height: float = 0.5,
        track_align_type="left",
        feature_track_ratio: float = 0.25,
        link_track_ratio: float = 1,
        theme="light",
        show_axis: bool = False,
        show_scale_bar: bool = False,
        show_scale_xticks: bool = False,
        font_size: int = 20,
        line_width: float = 1,
    ):
        """
        Parameters
        ----------
        seqs
            List of sequences
        seq_names, optional
            Names of the sequences to be overwitten, by default None
        """
        super().__init__(
            fig_width=fig_width,
            fig_track_height=fig_track_height,
            track_align_type=track_align_type,
            feature_track_ratio=feature_track_ratio,
            link_track_ratio=link_track_ratio,
            theme=theme,
            show_axis=show_axis,
        )

        if seq_names is None:
            seq_names = [seq.name for seq in seqs]

        for name, seq in zip(seq_names, seqs):
            self.add_feature_track(
                name,
                seq.length,
                labelsize=font_size,
                line_kws=dict(lw=line_width),
            )

        if show_scale_bar:
            self.set_scale_bar(labelsize=font_size)
        if show_scale_xticks:
            self.set_scale_xticks(labelsize=font_size)

        # Subtrack data
        self.subtrack_data = []

        # Marker data (for plotting a very limited number of points)
        self.marker_data = []

        # Line data (for plotting a very limited number of points)
        self.line_data = []

    def add_feature(
        self,
        bed_records_by_name: Dict[str, Sequence[BedRecord]],
        col: str,
        plotstyle="box",
        line_width=0.5,
        strand_plus_char: str = "+",
    ):
        """Add feature tracks (annotation on the sequences) to the plot.

        Parameters
        ----------
        bed_records_by_name
            Dict of {seq_name: [BedRecord]}.
            If a record `x` has `x.strand`, it will be considered.
        col
            Of the annotations
        plotstyle, optional
            By default, "box" if `self.box` is False, "bigbox" if `self.box` is True
        line_width, optional
            by default 0.5
        """
        for track in self.feature_tracks:
            if track.name not in bed_records_by_name:
                continue
            for x in bed_records_by_name[track.name]:
                track.add_feature(
                    x.b,
                    x.e,
                    1 if hasattr(x, "strand") and x.strand == strand_plus_char else -1,
                    plotstyle=plotstyle,
                    color=col,
                    lw=line_width,
                )

    def add_aln(self, df_by_pair: Dict[str, pd.DataFrame], curve=True, line_width=0.1):
        """Add alignment tracks to the plot.

        Parameters
        ----------
        bed_records_by_pair
            Dict of {(A_seq_name, B_seq_name): pd.DataFrame}}
        """
        for (a_name, b_name), df in df_by_pair.items():
            for _, row in df.iterrows():
                col = self.color_map.get(row["Type"], "black")
                is_inv = row["Type"] in ("INV", "INVTR", "INVDP")
                self.add_link(
                    (a_name, row["RefStart"], row["RefEnd"]),
                    (
                        b_name,
                        row["QryStart"] if not is_inv else row["QryEnd"],
                        row["QryEnd"] if not is_inv else row["QryStart"],
                    ),
                    color=col,
                    inverted_color=col,
                    lw=line_width,
                    curve=curve,
                )

    def add_subtrack(
        self,
        bed_records_by_name: Dict[str, Sequence[BedRecord]],
        x_func: Callable = lambda record: (record.b + record.e) // 2,
        y_func: Callable = lambda record: record.value,
        name=None,
        col="gray",
        ylim=(0, 100),
        ratio=1,
    ):
        """Add subtrack definition to the plot. The actual plot of the subtrack is
        not generated in this function but will be later generated in `show()`.

        Parameters
        ----------
        bed_records_by_name
            Dict of {seq_name: [BedRecord]}
        x_func, optional
            Mapping of BedRecord -> x-axis coordinate,
            by default lambda record: (record.b + record.e) // 2
        y_func, optional
            Mapping of BedRecord -> y-axis coordinate,
            by default lambda record: record.value
        """
        for track in self.feature_tracks:
            track.add_subtrack(name=name, ylim=ylim, ratio=ratio)
        self.subtrack_data.append((name, bed_records_by_name, x_func, y_func, col))

    def add_markers(
        self,
        bed_records_by_name: Dict[str, Sequence[BedRecord]],
        x_func: Callable = lambda record: (record.b + record.e) // 2,
        y_coords_by_name: Union[float, Dict[str, float]] = 1.5,
        col="gray",
        marker_size=15,
    ):
        """Add markers to the plot at a specific y-axis level for each seq.

        Parameters
        ----------
        bed_records_by_name
            Dict of {seq_name: [BedRecord]}
        x_func, optional
            Mapping of BedRecord -> x-axis coordinate,
            by default lambda record: (record.b + record.e) // 2
        y_coords_by_name, optional
            The level of y-axis the markers are plotted, by default 1.5
        """
        if isinstance(y_coords_by_name, float):
            y_coords_by_name = {
                track.name: y_coords_by_name for track in self.feature_tracks
            }

        self.marker_data.append(
            (bed_records_by_name, x_func, y_coords_by_name, col, marker_size)
        )

    def add_lines(
        self,
        bed_records_by_name: Dict[str, Sequence[BedRecord]],
        y_coords_by_name: Union[float, Dict[str, float]] = 1.5,
        col="gray",
        line_width=1,
    ):
        """Add markers to the plot at a specific y-axis level for each seq.

        Parameters
        ----------
        bed_records_by_name
            Dict of {seq_name: [BedRecord]}
        x_func, optional
            Mapping of BedRecord -> x-axis coordinate,
            by default lambda record: (record.b + record.e) // 2
        y_coords_by_name, optional
            The level of y-axis the markers are plotted, by default 1.5
        """
        if isinstance(y_coords_by_name, float):
            y_coords_by_name = {
                track.name: y_coords_by_name for track in self.feature_tracks
            }

        self.line_data.append((bed_records_by_name, y_coords_by_name, col, line_width))

    def show(self, dpi=300, out_image=None):
        """Show the plot.

        Parameters
        ----------
        out_image, optional
            Path of the out image file, by default None
            The format is e.g.:
              - `out.svg` for single output
              - `out.{svg,pdf,png}` for multiple outputs
        """
        # Subtrack axes are generated by this
        fig = self.plotfig(dpi=dpi)

        # Plot the subtrack data
        for track in self.feature_tracks:
            for name, bed_records_by_name, x_func, y_func, col in self.subtrack_data:
                if track.name not in bed_records_by_name:
                    continue
                data = bed_records_by_name[track.name]
                x, y = list(map(x_func, data)), list(map(y_func, data))
                subtrack = track.get_subtrack(name)
                subtrack.ax.fill_between(subtrack.transform_coord(x), y, color=col)

        # Plot the marker data
        for track in self.feature_tracks:
            for (
                bed_records_by_name,
                x_func,
                y_coords_by_name,
                col,
                marker_size,
            ) in self.marker_data:
                if track.name not in bed_records_by_name:
                    continue
                data = bed_records_by_name[track.name]
                x = list(map(x_func, data))
                track.ax.scatter(
                    x=track.transform_coord(x),
                    y=[y_coords_by_name[track.name]] * len(x),
                    color=col,
                    s=marker_size,
                    clip_on=False,
                )

        # Plot the line data
        for track in self.feature_tracks:
            for (
                bed_records_by_name,
                y_coords_by_name,
                col,
                line_width,
            ) in self.line_data:
                if track.name not in bed_records_by_name:
                    continue
                data = bed_records_by_name[track.name]
                x_b, x_e = [x.b for x in data], [x.e for x in data]
                track.ax.hlines(
                    xmin=track.transform_coord(x_b),
                    xmax=track.transform_coord(x_e),
                    y=[y_coords_by_name[track.name]] * len(data),
                    color=col,
                    lw=line_width,
                    clip_on=False,
                )

        plt.figure(fig)
        plt.show(fig)

        if out_image is not None:
            if out_image.endswith("}"):  # multiple output formats
                data = out_image[:-1].split("{")
                assert len(data) == 2, f"Invalid format: {out_image}"
                prefix = data[0]
                exts = data[1].split(",")
                out_fnames = [f"{prefix}{ext}" for ext in exts]
            else:  # single output formats
                out_fnames = [out_image]

            for out_fname in out_fnames:
                fig.savefig(out_fname)
