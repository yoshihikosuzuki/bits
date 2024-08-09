from typing import List, Optional, Sequence, Tuple, Union

import plotly.graph_objects as go
import plotly_light as pl
import pysam

from .._io import load_bam
from .._type import BedRecord, SegRecord


def align_pileup(
    read_start_ends: List[Tuple[int, int]], min_spacing: int = 100
) -> List[int]:
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
            if max_e + min_spacing < b:
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
    in_bam_or_reads: Union[str, Sequence[pysam.AlignedSegment]],
    region: Union[str, SegRecord],
    require_span: bool = False,
    min_read_len: Optional[int] = None,
    min_spacing: Union[float, int] = 0.01,
    show_suppl: bool = False,
    color_strand: bool = False,
    marker_size: int = 2,
    line_width: int = 2,
    line_width_clip: int = 1,
    width: int = 1000,
    height: int = 400,
    layout: Optional[go.Layout] = None,
    use_webgl: bool = False,
    return_fig: bool = False,
) -> Optional[go.Figure]:
    """Show the read pileup plot given a .bam file.

    Parameters
    ----------
    in_bam_or_reads
        Path to a .bam file of mapped reads. Or, a list of mappings.
    region
        A string (e.g. "chr1:1000-2000") or BedRecord.
    require_span
        Show only reads spanning `region`.
    min_read_len
        To be shown
    min_spacing
        Min size of gaps between read bars. Ratio w.r.t. `region` size (in float, between 0 and 1)
        or base pairs (in int)
    show_suppl
        Show supplementary alignments and secondary alignments as well as primary alignments.
    color_strand
        Use different colors for different strands.
    return_fig, optional
        Return pl.Figure object instead of drawing a plot, by default False

    Returns
    -------
        pl.Figure if `return_fig` is True, otherwise None
    """
    region = SegRecord.from_string(region) if isinstance(region, str) else region
    min_spacing = (
        region.length * min_spacing
        if 0 < min_spacing and min_spacing < 1
        else min_spacing
    )

    ##### Preparing necassary data #####
    if isinstance(in_bam_or_reads, str):
        reads = load_bam(in_bam_or_reads, region, require_span)
    else:
        reads = in_bam_or_reads

    if not show_suppl:
        reads = list(filter(lambda read: read.flag in (0, 16), reads))
    if min_read_len is not None:
        reads = list(filter(lambda read: read.query_length >= min_read_len, reads))

    # soft-clipped length for each end
    clip_lens = [(read.qstart, read.inferred_length - read.qend) for read in reads]

    # calculate row ID for each read
    read_rows = align_pileup(
        [
            (read.reference_start - b_clip, read.reference_end + e_clip)
            for read, (b_clip, e_clip) in zip(reads, clip_lens)
        ],
        min_spacing=min_spacing,
    )

    ##### Filter functions for reads #####
    def _filter_by_read(cond_read):
        return list(
            zip(*(filter(lambda t: cond_read(t[0]), zip(reads, read_rows, clip_lens))))
        )

    def _filter_by_flag(
        is_primary: Optional[bool] = None, is_forward: Optional[bool] = None
    ):
        assert is_primary in (None, True, False)
        assert is_forward in (None, True, False)
        if is_primary is None and is_forward is None:
            return reads, read_rows, clip_lens
        return _filter_by_read(
            lambda read: (
                True
                if is_primary is None
                else read.flag & 2048 == (0 if is_primary else 2048)
            )
            and (
                True
                if is_forward is None
                else read.flag & 16 == (0 if is_forward else 16)
            )
        )

    ##### Defining colors #####
    ## Colors for aligned segments and annotations
    # {(is_primary, is_forward): color}
    COL_TABLE = {
        (True, True): pl.colors["blue"],
        (True, False): pl.colors["yellow"],
        (False, True): pl.colors["darkblue"],
        (False, False): pl.colors["darkyellow"],
        (None, True): pl.colors["blue"],
        (None, False): pl.colors["yellow"],
        (True, None): "gray",
        (False, None): "gray",
        (None, None): "gray",
    }
    ## Colors for soft_clipped regions
    # {(is_primary: color}
    COL_TABLE_READ = {
        True: pl.colors["red"],
        False: pl.colors["darkred"],
        None: pl.colors["red"],
    }

    ##### Making traces #####
    ## Hover text
    def _trace_text(
        is_primary: Optional[bool], is_forward: Optional[bool], which_end: str
    ):
        _reads, _read_rows, _clip_lens = _filter_by_flag(is_primary, is_forward)
        return pl.scatter(
            x=[
                read.reference_start if which_end == "start" else read.reference_end
                for read in _reads
            ],
            y=[row_id for row_id in _read_rows],
            text=[
                "<br>".join(
                    [
                        f"{read.qname}",
                        f"flag = {read.flag}",
                        f"ref: {read.reference_start:8,} - {read.reference_end:8,}",
                        f"read: {read.qstart:6,} - {read.qend:6,}",
                        f"clip: {b_clip:6,} bp and {e_clip:6,} bp",
                    ]
                )
                for read, (b_clip, e_clip) in zip(_reads, _clip_lens)
            ],
            col=COL_TABLE[(is_primary, is_forward)],
            marker_size=marker_size,
            use_webgl=use_webgl,
        )

    traces_text = [
        _trace_text(
            is_primary,
            is_forward,
            which_end,
        )
        for is_primary in ((True, False) if show_suppl else (None,))
        for is_forward in ((True, False) if color_strand else (None,))
        for which_end in ("start", "end")
    ]

    ## Entire read
    def _trace_read(is_primary: Optional[bool]):
        _reads, _read_rows, _clip_lens = _filter_by_flag(is_primary=is_primary)
        return pl.lines(
            [
                (
                    read.reference_start - b_clip,
                    row_id,
                    read.reference_end + e_clip,
                    row_id,
                )
                for read, row_id, (b_clip, e_clip) in zip(
                    _reads, _read_rows, _clip_lens
                )
            ],
            col=COL_TABLE_READ[is_primary],
            width=line_width_clip,
            use_webgl=use_webgl,
        )

    traces_read = [
        _trace_read(is_primary)
        for is_primary in ((True, False) if show_suppl else (None,))
    ]

    ## Aligned segment
    def _trace_align(is_primary: Optional[bool], is_forward: Optional[bool]):
        _reads, _read_rows, _ = _filter_by_flag(
            is_primary=is_primary, is_forward=is_forward
        )
        return pl.lines(
            [
                (read.reference_start, row_id, read.reference_end, row_id)
                for read, row_id in zip(_reads, _read_rows)
            ],
            col=COL_TABLE[(is_primary, is_forward)],
            width=line_width,
            use_webgl=use_webgl,
        )

    traces_align = [
        _trace_align(
            is_primary,
            is_forward,
        )
        for is_primary in ((True, False) if show_suppl else (None,))
        for is_forward in ((True, False) if color_strand else (None,))
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
                title=f"{region.to_string(comma=True)}",
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
