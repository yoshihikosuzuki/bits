from typing import List, Tuple

import plotly_light as pl
import pysam


def align_pileup(read_start_ends: List[Tuple[int, int]], margin=100):
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


def show_read_pileup(in_bam: str, region: str, return_fig: bool = False):
    chrom, b_e = region.split(":")
    b, e = b_e.split("-")
    b, e = int(b), int(e)

    reads = [read for read in pysam.AlignmentFile(in_bam, "rb").fetch(chrom, b, e)]
    read_rows = align_pileup(
        [(read.reference_start, read.reference_end) for read in reads]
    )
    traces = pl.lines(
        [
            (read.reference_start, row_id, read.reference_end, row_id)
            for read, row_id in zip(reads, read_rows)
        ],
        col="blue",
    )
    layout = pl.layout(
        width=1300,
        height=400,
        x_bounding_line=True,
        y_bounding_line=True,
        x_mirror=True,
        y_mirror=True,
        x_grid=True,
        y_grid=True,
        grid_col="#eeeeee",
        x_title=f"Genomic Position on {chrom}",
        y_title="Reads",
        x_ticks_minor="outside",
        x_tickformat=",.0f",
    )
    fig = pl.figure(traces, layout)
    if return_fig:
        return fig
    pl.show(fig)
