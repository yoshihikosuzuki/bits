# NOTE: Visualization for whole genome annotation data

import plotly_light as pl


def gen_chr_plot(refs, names=None, title=None, width=700, height=1000):
    if names is None:
        names = [ref.name for ref in refs]
    return pl.figure(
        pl.scatter(
            [0 for ref in refs for x in ["t", "b"]],
            [f"{ref.name}_{x}" for ref in refs for x in ["t", "b"]],
            text=[
                name if x == "t" else ""
                for ref, name in zip(refs, names)
                for x in ["t", "b"]
            ],
            mode="text",
            text_pos="top right",
            text_size=12,
            use_webgl=False,
        ),
        pl.layout(
            width=width,
            height=height,
            title=title,
            x_zeroline=False,
            y_category=True,
            y_reversed=True,
            y_ticklabel=False,
            shapes=[
                pl.rect(
                    0,
                    f"{ref.name}_b",
                    ref.length,
                    f"{ref.name}_t",
                    layer="below",
                    fill_col="white",
                    frame_width=0.5,
                )
                for ref in refs
            ],
        ),
    )


def add_bed_records(fig, refs, bed_records, col, opacity=1, width=1):
    names = set([ref.name for ref in refs])
    bed_records = list(filter(lambda x: x.chr in names, bed_records))
    fig.add_trace(
        pl.scatter(
            [x.b for x in bed_records],
            [f"{x.chr}_t" for x in bed_records],
            text=[f"{x.chr}:{x.b:,}-{x.e:,}<br>{x.e - x.b:,} bp" for x in bed_records],
            col=col,
            opacity=0,
        )
    )
    fig.add_trace(
        pl.scatter(
            [x.e for x in bed_records],
            [f"{x.chr}_b" for x in bed_records],
            text=[f"{x.chr}:{x.b:,}-{x.e:,}<br>{x.e - x.b:,} bp" for x in bed_records],
            col=col,
            opacity=0,
        )
    )
    shapes = [
        pl.rect(
            x.b,
            f"{x.chr}_b",
            x.e,
            f"{x.chr}_t",
            layer="above",
            opacity=opacity,
            fill_col=col,
            frame_col=col,
            frame_width=width,
        )
        for x in bed_records
    ]
    fig.layout.shapes = tuple(list(fig.layout.shapes) + shapes)
