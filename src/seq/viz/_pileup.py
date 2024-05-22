from typing import List, Optional, Tuple

import plotly.graph_objects as go
import plotly_light as pl

from .._io import load_bam


def align_pileup(read_start_ends: List[Tuple[int, int]], margin=100) -> List[int]:
    """Given a list of (start, end) intervals, calculate the min row ID for each
    interval s.t. each interval does not overlap with other intervals.

    Parameters
    ----------
    read_start_ends
        A list of `(interval_start, interval_end)`.
    margin, optional
        Adjacent intervals must be separated by this value, by default 100

    Returns
    -------
        A list of row ID for each interval.
    """
    read_rows = [None] * len(read_start_ends)
    row_max_es = []
    for read_id, (b, e) in enumerate(read_start_ends):
        read_row = None
        for row_id, max_e in enumerate(row_max_es):
            if max_e + margin < b:
                read_row = row_id
                break
        if read_row is None:
            read_row = len(row_max_es)
            row_max_es.append(e)
        else:
            row_max_es[read_row] = e
        read_rows[read_id] = read_row
    return read_rows


def show_read_pileup(
    in_bam: str,
    region: str,
    show_suppl: bool = False,
    color_strand: bool = False,
    color_suppl: bool = False,
    marker_size: int = 2,
    line_width: int = 2,
    width: int = 1000,
    height: int = 400,
    layout: Optional[go.Layout] = None,
    return_fig: bool = False,
) -> Optional[go.Figure]:
    """Show the read pileup plot given a .bam file.

    Parameters
    ----------
    in_bam
        Path to a .bam file of mapped reads.
    region
        e.g. "chr1:1000-2000"
    show_suppl
        Show supplementary alignments and secondary alignments as well as primary alignments.
    color_strand
        Use different colors for different strands.
    color_suppl
        Use different colors for primary alignments and other alignments.
    return_fig, optional
        Return pl.Figure object instead of drawing a plot, by default False

    Returns
    -------
        pl.Figure if `return_fig` is True, otherwise None
    """
    reads = load_bam(in_bam, region)
    if not show_suppl:
        reads = list(filter(lambda read: read.flag in (0, 16), reads))

    # soft-clipped length for each end
    sc_lens = [(read.qstart, read.inferred_length - read.qend) for read in reads]

    # calculate row ID for each read
    read_rows = align_pileup(
        [
            (read.reference_start - b_sc, read.reference_end + e_sc)
            for read, (b_sc, e_sc) in zip(reads, sc_lens)
        ]
    )

    # Hover text
    traces_text = [
        pl.scatter(
            x=x,
            y=[row_id for row_id in read_rows],
            text=[
                "<br>".join(
                    [
                        f"{read.qname}",
                        f"flag = {read.flag}",
                        f"ref: {read.reference_start:8,} - {read.reference_end:8,}",
                        f"read: {read.qstart:6,} - {read.qend:6,}",
                        f"clip: {b_sc:6,} bp and {e_sc:6,} bp",
                    ]
                )
                for read, (b_sc, e_sc) in zip(reads, sc_lens)
            ],
            col=(
                "gray"
                if not color_strand
                else [
                    pl.colors["blue"] if read.flag == 0 else pl.colors["yellow"]
                    for read in reads
                ]
            ),
            marker_size=marker_size,
        )
        for x in (
            [read.reference_start for read in reads],  # 5'-end
            [read.reference_end for read in reads],  # 3'-end
        )
    ]
    # Entire read
    traces_read = [
        pl.lines(
            [
                (
                    read.reference_start - b_sc,
                    row_id,
                    read.reference_end + e_sc,
                    row_id,
                )
                for read, row_id, (b_sc, e_sc) in zip(reads, read_rows, sc_lens)
            ],
            col=pl.colors["red"],
            width=line_width,
        )
    ]
    # Aligned segment
    traces_align = [
        pl.lines(
            [
                (read.reference_start, row_id, read.reference_end, row_id)
                for read, row_id in zip(reads, read_rows)
                if not color_strand or read.flag == flag
            ],
            col=col,
            width=line_width,
        )
        for flag, col in (
            ((0, pl.colors["blue"]), (16, pl.colors["yellow"]))
            if color_strand
            else ((None, "gray"),)
        )
    ]

    fig = pl.figure(
        traces_text + traces_read + traces_align,
        pl.merge_layout(
            pl.layout(
                width=width,
                height=height,
                x_bounding_line=True,
                y_bounding_line=True,
                x_mirror=True,
                y_mirror=True,
                x_grid=True,
                y_grid=True,
                grid_col="#eeeeee",
                title=f"{region}",
                x_title="Genomic Position",
                y_title="Reads",
                x_ticks_minor="outside",
                x_tickformat=",.0f",
            ),
            layout,
        ),
    )
    if return_fig:
        return fig
    pl.show(fig)
