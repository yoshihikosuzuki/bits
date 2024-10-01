from dataclasses import dataclass
from typing import Optional, Sequence, Union

import plotly.graph_objects as go
import plotly_light as pl

from .._io import load_fasta
from .._type import BedRecord, FastaRecord


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
