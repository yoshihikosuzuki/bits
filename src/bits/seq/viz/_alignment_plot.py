from typing import Optional

import bits.util as bu
import plotly.graph_objects as go
import plotly_light as pl


def show_alignment_plot(
    in_paf: str, width: int = 600, height: int = 600, layout: Optional[go.Layout] = None
):
    """Show an alignment plot.
    Currently there are many assumptions:
        - only single query vs single target is supported
        - only forward vs forward alignment is supported
        - aspect ratio is not conserved

    Parameters
    ----------
    in_paf
        Input .paf file between two sequences.
    """
    data = list(
        map(
            lambda line: line.split("\t"),
            bu.run_command(f"cat {in_paf}").strip().split("\n"),
        )
    )

    ql = int(data[0][1])
    tl = int(data[0][6])

    pl.show(
        pl.lines(
            [(int(x[2]), int(x[7]), int(x[3]), int(x[8])) for x in data],
        ),
        pl.layout(
            width=width + 50,
            height=height / int(ql) * int(tl) + 50,
            x_range=(0, ql),
            y_range=(0, tl),
            x_zeroline=True,
            y_zeroline=True,
            x_title="Query",
            y_title="Target",
            x_grid=False,
            y_grid=False,
            x_bounding_line=True,
            x_mirror=True,
            y_bounding_line=True,
            y_mirror=True,
            anchor_axes=True,
        ),
    )
