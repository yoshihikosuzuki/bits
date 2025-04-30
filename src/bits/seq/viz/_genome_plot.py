from typing import Callable, Dict, Sequence

import matplotlib.pyplot as plt
from pygenomeviz import GenomeViz

from .._type import BedRecord, SeqRecord


class GenomePlotGV(GenomeViz):
    """Wrapper of pyGenomeViz for genome plot.

    Example usage:
    ```
    > seqs = bs.load_fasta("in.fasta")
    > anno = bs.load_bed("anno.bed", by_chrom=True)
    > cov = bs.load_bed("cov.bedgraph", attrs=[("cov", float)], by_chrom=True)
    >
    > gv = bs.GenomePlotGV(seqs, show_scale_xticks=True, box=True)
    > gv.add_feature(anno)
    > gv.add_subtrack(cov)
    > gv.show()
    ```
    """

    def __init__(
        self,
        seqs: Sequence[SeqRecord],
        fig_width: float = 5,
        fig_track_height: float = 0.1,
        track_align_type="left",
        feature_track_ratio: float = 2,
        link_track_ratio: float = 1,
        theme="light",
        show_axis: bool = False,
        show_scale_bar: bool = False,
        show_scale_xticks: bool = False,
        font_size: int = 10,
        box: bool = False,
        line_width: float = 0.5,
    ):
        """_summary_

        Parameters
        ----------
        seqs
            A list of sequences to be plotted.
        show_scale_bar, optional
            If True, draw a scale bar at the bottom right, by default False
        show_scale_xticks, optional
            If True, draw bounding line, ticks, and tick labels on x-axis , by default False
        font_size, optional
            Font size of sequence names and axis labels, by default 10
        box, optional
            If True, draw each sequence as a box instead of bar, by default False
        line_width, optional
            Of boxes (valid only when `box` is True), by default 0.5

        (Others: same as pygenomeviz.GenomeViz)
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
        self.box = box

        if show_scale_bar:
            self.set_scale_bar(labelsize=font_size)
        if show_scale_xticks:
            self.set_scale_xticks(labelsize=font_size)

        # Sequence tracks definition
        for seq in seqs:
            track = self.add_feature_track(
                name=seq.name,
                segments=seq.length,
                labelsize=font_size,
                line_kws=dict(lw=0) if box is True else None,
            )
            if self.box is True:
                track.add_feature(
                    0,
                    seq.length,
                    plotstyle="bigbox",
                    fc="none",
                    ec="black",
                    lw=line_width,
                )

        # Subtrack data
        self.subtrack_data = []

    def add_feature(
        self,
        bed_records_by_name: Dict[str, Sequence[BedRecord]],
        col: str,
        plotstyle=None,
        line_width=0.5,
        strand_plus_char: str = "+",
    ):
        """Add feature tracks (annotation on the sequences) to the plot.

        Parameters
        ----------
        bed_records_by_name
            Dict of {seq_name: [BedRecord]}
        col
            Of the annotations
        plotstyle, optional
            By default, "box" if `self.box` is False, "bigbox" if `self.box` is True
        line_width, optional
            by default 0.5
        """
        if plotstyle is None:
            plotstyle = "bigbox" if self.box is True else "box"
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
