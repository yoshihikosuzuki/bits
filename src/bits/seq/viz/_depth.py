from typing import Optional, Sequence

import plotly.graph_objects as go
import plotly_light as pl

from .._type import BedRecord
from ._region_feature import trace_bed_attr


def fig_depth_detail(
    data: Sequence[BedRecord],
    mean_cov: Optional[float] = None,
    line_width: float = 2,
    layout: Optional[go.Layout] = None,
) -> go.Figure:
    """Coverage data with detailed information about min/max/median coverage to a Figure object.

    @ data
        A list of coverage records that have `.b`, `.min`, `.max`, `.med` as variables.
    @ mean_cov
        Global mean coverage.
    """
    return pl.figure(
        (
            [
                # global mean coverage
                pl.lines_shape(
                    (0, mean_cov, 1, mean_cov),
                    width=line_width,
                    col="gray",
                    name="global mean",
                )
            ]
            if mean_cov is not None
            else []
        )
        + [
            # min coverage per bin
            trace_bed_attr(
                data, "min", pl.colors["red"], line_width=line_width, show_legend=True
            ),
            # max coverage per bin
            trace_bed_attr(
                data,
                "max",
                pl.colors["yellow"],
                line_width=line_width,
                show_legend=True,
            ),
            # median coverage per bin
            trace_bed_attr(
                data,
                "med",
                pl.colors["blue"],
                line_width=line_width,
                name="median",
                show_legend=True,
            ),
        ],
        pl.merge_layout(
            pl.layout(
                x_bounding_line=True,
                y_bounding_line=True,
                x_mirror=True,
                y_mirror=True,
                x_grid=True,
                y_grid=True,
            ),
            layout,
        ),
    )


def fig_depth_mean(
    data: Sequence[BedRecord],
    mean_cov: Optional[float] = None,
    line_width: float = 1,
    col: str = "gray",
    layout: Optional[go.Layout] = None,
) -> go.Figure:
    """Data of some coverage value per bin to a Figure object.

    @ data
        A list of coverage records that have `.b`, `.cov` as variables.
    @ mean_cov
        Global mean coverage.
    """
    return pl.figure(
        (
            [
                # global mean coverage
                pl.lines_shape(
                    (0, mean_cov, 1, mean_cov),
                    width=line_width,
                    col="gray",
                    name="global mean",
                )
            ]
            if mean_cov is not None
            else []
        )
        + [trace_bed_attr(data, attr="cov", col=col, line_width=line_width)],
        pl.merge_layout(
            pl.layout(
                x_bounding_line=True,
                y_bounding_line=True,
                x_mirror=True,
                y_mirror=True,
                x_grid=True,
                y_grid=True,
            ),
            layout,
        ),
    )
