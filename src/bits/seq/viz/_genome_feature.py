## NOTE: deprecated. use GenomePlotGV instead.

from dataclasses import dataclass
from typing import Optional, Sequence, Union

import plotly.graph_objects as go
import plotly_light as pl

from .._io import load_fasta
from .._type import BedRecord, FastaRecord

## NOTE: seq に番号をつけた版

# def gen_chr_plot(seqs, names=None, title=None, width=750, height=1000):
#     if names is None:
#         names = [seq.name for seq in seqs]

#     seq_to_y = {seq.name: i for i, seq in enumerate(reversed(seqs), start=1)}

#     return pl.figure(
#         pl.rects([
#             (
#                 0,
#                 seq_to_y[seq.name] - 0.25,
#                 seq.length,
#                 seq_to_y[seq.name] + 0.25
#             ) for seq in seqs],
#             fill_col=None,
#         ),
#         pl.layout(
#             width=width,
#             height=height,
#             title=title,
#             box=False,
#             x_bounding_line=True,
#             x_ticks=True
#         ).update(dict(yaxis=dict(
#             tickvals=list(seq_to_y.values()),
#             ticktext=list(seq_to_y.keys()),
#             range=[min(seq_to_y.values()) - 0.5, max(seq_to_y.values()) + 0.5],
#         )))
#     )


# def add_bed_records(fig, seqs, bed_records, col, name="", show_legend=False, opacity=1, width=1):
#     seq_to_y = {seq.name: i for i, seq in enumerate(reversed(seqs), start=1)}
#     names = set([seq.name for seq in seqs])
#     bed_records = list(filter(lambda x: x.chr in names, bed_records))

#     fig.add_trace(
#         pl.rects(
#             [(x.b, seq_to_y[x.chr] - 0.25, x.e, seq_to_y[x.chr] + 0.25)
#              for x in bed_records],
#             opacity=opacity,
#             fill_col=col,
#             frame_col=col,
#             frame_width=width,
#             name=name,
#         )
#     )


## NOTE: 旧版


@dataclass
class GenomePlot:
    """Class for a genome-wide plot."""

    seqs: Union[str, Sequence[FastaRecord]]
    names: Optional[Sequence[str]] = None
    layout: Optional[go.Layout] = None

    def __post_init__(self):
        if isinstance(self.seqs, str):
            self.seqs = load_fasta(self.seqs)
        if self.names is None:
            self.names = [seq.name for seq in self.seqs]

        self.fig = pl.figure(
            pl.scatter(
                [0 for seq in self.seqs for x in ["t", "b"]],
                [f"{seq.name}_{x}" for seq in self.seqs for x in ["t", "b"]],
                text=[
                    name if x == "t" else ""
                    for seq, name in zip(self.seqs, self.names)
                    for x in ["t", "b"]
                ],
                mode="text",
                text_pos="top right",
                text_size=12,
                use_webgl=False,
            ),
            pl.merge_layout(
                pl.layout(
                    width=700,
                    height=1000,
                    x_zeroline=False,
                    y_category=True,
                    y_reversed=True,
                    y_ticklabel=False,
                    shapes=[
                        pl.rect(
                            0,
                            f"{seq.name}_b",
                            seq.length,
                            f"{seq.name}_t",
                            layer="below",
                            fill_col="white",
                            frame_width=0.5,
                        )
                        for seq in self.seqs
                    ],
                ),
                self.layout,
            ),
        )

    def add_bed_records(
        self,
        bed_records: Sequence[BedRecord],
        width: float = 1,
        col: str = "black",
        opacity: Optional[float] = None,
        layer: str = "above",
        name: Optional[str] = None,
        show_legend: bool = False,
        use_webgl: bool = False,
    ):
        names_set = set(self.names)
        bed_records = list(filter(lambda x: x.chr in names_set, bed_records))

        # TODO: use scatter plot with some line_width instead of shapes + dummpy scatter???

        # for hover text at start positions
        self.fig.add_trace(
            pl.scatter(
                [x.b for x in bed_records],
                [f"{x.chr}_t" for x in bed_records],
                text=[
                    f"{x.chr}:{x.b:,}-{x.e:,}<br>{x.e - x.b:,} bp" for x in bed_records
                ],
                col=col,
                opacity=0,
            )
        )
        # for hover text at end positions
        self.fig.add_trace(
            pl.scatter(
                [x.e for x in bed_records],
                [f"{x.chr}_b" for x in bed_records],
                text=[
                    f"{x.chr}:{x.b:,}-{x.e:,}<br>{x.e - x.b:,} bp" for x in bed_records
                ],
                col=col,
                opacity=0,
            )
        )
        # for legend
        self.fig.add_trace(
            pl.scatter(
                [None],
                [None],
                marker_size=10,
                col=col,
                opacity=1,
                name=name,
                show_legend=show_legend,
                use_webgl=use_webgl,
            )
        )
        self.fig.layout.shapes = tuple(
            list(self.fig.layout.shapes)
            + [
                pl.rect(
                    x.b,
                    f"{x.chr}_b",
                    x.e,
                    f"{x.chr}_t",
                    layer=layer,
                    opacity=1,
                    fill_col=col,
                    frame_col=col,
                    frame_width=width,
                )
                for x in bed_records
            ]
        )
        return self

    def show(self) -> Optional[go.Figure]:
        pl.show(self.fig)
