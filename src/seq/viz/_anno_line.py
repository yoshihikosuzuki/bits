from typing import Optional, Sequence, Union

import plotly.graph_objects as go
import plotly_light as pl

from .._io import load_bed
from .._type import BedRecord


def trace_depth_detail(
    data: Sequence[BedRecord],
    mean_cov: float,
    max_cov: float,
    chrom_len: int,
    line_width: float = 2,
    title: Optional[str] = None,
) -> go.Figure:
    """Coverage data with detailed information about min/max/median coverage to a Figure object.

    Parameters
    ----------
    data
        A list of coverage records that have `.b`, `.min`, `.max`, `.med` as variables.
    mean_cov
        Global mean coverage.
    max_cov
        Max coverage to be shown.
    chrom_len
        Sequence length.
    line_width, optional
        by default 2
    title, optional
        by default None

    Returns
    -------
        Figure object
    """
    x = [r.b for r in data]
    return pl.figure(
        [
            pl.lines(
                (0, mean_cov, chrom_len, mean_cov),
                width=line_width,
                col="gray",
                name="global mean",
                use_webgl=False,
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.min for r in data],
                mode="lines",
                col=pl.colors["red"],
                name="min",
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.max for r in data],
                mode="lines",
                col=pl.colors["yellow"],
                name="max",
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.med for r in data],
                mode="lines",
                col=pl.colors["blue"],
                name="median",
                line_width=line_width,
                show_legend=True,
            ),
        ],
        pl.layout(
            title=title,
            y_range=(0, max_cov),
            x_bounding_line=True,
            y_bounding_line=True,
            x_mirror=True,
            y_mirror=True,
            x_grid=True,
            y_grid=True,
        ),
    )


def trace_depth_mean(
    data: Sequence[BedRecord],
    mean_cov: float,
    max_cov: float,
    chrom_len: int,
    line_width: float = 1,
    col: str = "gray",
    title: Optional[str] = None,
) -> go.Figure:
    """_summary_

    Parameters
    ----------
    data
       Each record has to have `.b` and `.cov` variables.
    mean_cov
        _description_
    max_cov
        _description_
    chrom_len
        _description_
    title
        _description_

    Returns
    -------
        _description_
    """
    return pl.figure(
        [
            pl.lines((0, mean_cov, chrom_len, mean_cov), width=line_width, col=col),
            pl.scatter(
                [r.b for r in data],
                [r.cov for r in data],
                mode="lines",
                col=col,
                line_width=line_width,
                use_webgl=False,
            ),
        ],
        pl.layout(
            title=title,
            y_range=(0, max_cov),
            x_bounding_line=True,
            y_bounding_line=True,
            x_mirror=True,
            y_mirror=True,
            x_grid=True,
            y_grid=True,
        ),
    )


def trace_bed(
    data: Union[str, Sequence[BedRecord]],
    chrom: str,
    name: str,
    col: str,
    width: float = 8,
    pad: float = 0,
    return_trace: bool = True,
) -> Union[go.Figure, go.Trace]:
    """Make a Figure object from a .bed file or a list of BedRecords

    Parameters
    ----------
    data
        Name of a .bed file, or a list of `BedRecords`
    pad
        Size of paddings for each record. [`r.b - pad`..`r.b + pad`] is drawn.
    """
    if isinstance(data, str):
        data = load_bed(data)
    data = list(filter(lambda x: x.chr == chrom, data))

    if len(data) == 0:
        return pl.lines(
            [(0, name, 0, name)], width=0, col=col, name=name, show_legend=False
        )

    return pl.lines(
        [(max(x.b - pad, 0), name, x.e + pad, name) for x in data],
        text=[f"{x.chr}:{x.b}-{x.e}" for x in data],
        width=width,
        col=col,
        name=name,
        show_legend=False,
    )
