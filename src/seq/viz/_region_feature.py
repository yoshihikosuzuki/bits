from typing import Optional, Sequence, Type, Union

import plotly.graph_objects as go
import plotly_light as pl
from logzero import logger
from pysam import VariantRecord

from .._io import filter_bed, load_bed, load_vcf
from .._type import BedRecord, SegRecord


def trace_bed(
    data: Union[str, BedRecord, Sequence[BedRecord]],
    region: Optional[Union[str, SegRecord]] = None,
    text: Optional[Sequence[str]] = None,
    name: str = "",
    col: Optional[str] = None,
    width: float = 8,
    pad: float = 0,
    min_len: float = 0,
    use_webgl: bool = False,
) -> go.Trace:
    """For plotting the intervals of the bed records.

    NOTE: Add an empty `trace_bed` both at the beginning and at the enf of the list of
          `trace_bed` for spacing:

          > pl.show(
                [bs.trace_bed([], name=" ")] +      # dummy
                [bs.trace_bed(...) for ...] +       # main
                [be.trace_bed([], name="  ")],      # dummy
                layout=pl.layout(y_reversed=True)   # so track order becomes up to down
            )

    Parameters
    ----------
    data
        Name of a .bed file, or a list of `BedRecords`.
    region
        Chromosome name (str) or a region to be shown.
        Note that only single sequence is supported for this function.
    text
        Hover texts.
    text_above
        Texts shown on each record.
    name
        Track name.
    pad
        Size of paddings for each record. [`r.b - pad`..`r.b + pad`] is drawn.
    min_len
        Minimum length in the plot for each record.
    """
    if isinstance(data, str):
        data = load_bed(data, region)
    else:
        if isinstance(data, BedRecord):
            data = [data]
        data = filter_bed(data, region)

    if len(data) == 0:
        # Dummy trace just to let the track exist
        return pl.lines(
            [(0, name, 0, name)], width=0, col=col, name=name, show_legend=False
        )
    else:
        pads = [max((min_len - x.length) / 2, pad) for x in data]
        return pl.lines(
            [(max(x.b - p, 0), name, x.e + p, name) for x, p in zip(data, pads)],
            text=[f"{x.chr}:{x.b}-{x.e}" for x in data] if text is None else text,
            width=width,
            col=col,
            name=name,
            show_legend=False,
            use_webgl=use_webgl,
        )


def trace_bed_attr(
    data: Union[str, Sequence[BedRecord]],
    attr: str,
    attr_type: Optional[Type] = None,
    region: Optional[Union[str, SegRecord]] = None,
    col: Optional[str] = None,
    line_width: float = 2,
    name: Optional[str] = None,
    show_legend: bool = False,
    use_webgl: bool = False,
) -> go.Trace:
    """For plotting an attribute of the bed(graph) records.

    Parameters
    ----------
    data
        Name of a .bed(graph) file, or a list of `BedRecords`.
    attr
        Name of the attribute to be shown in the y-axis.
    attr_type
        Type of the attribute to be shown in the y-axis. Necessary if `data` is str.
    region
        Chromosome name (str) or a region to be shown.
        Note that only single sequence is supported for this function.
    name
        Track name.
    """
    if isinstance(data, str):
        assert (
            attr_type is not None
        ), "`attr_type` needs to be specified when `data` is str."
        data = load_bed(data, region, attrs=[(attr, attr_type)])
    else:
        if attr_type is not None:
            logger.info("Ignored `attr_type` because pre-loaded data is specified.")
        if isinstance(data, BedRecord):
            data = [data]
        data = filter_bed(data, region)

    if name is None:
        name = attr
    return pl.scatter(
        [r.b for r in data],
        [getattr(r, attr) for r in data],
        mode="lines",
        col=col,
        name=name,
        line_width=line_width,
        show_legend=show_legend,
        use_webgl=use_webgl,
    )


def trace_vcf(
    data: Union[str, Sequence[VariantRecord]],
    chrom: Optional[str] = None,
    name: str = "track",
    col: Optional[str] = None,
    width: float = 8,
    pad: float = 0,
    single_sample: bool = False,
    use_webgl: bool = False,
):
    """Make a Figure object from a .vcf file or a list of vcf.model._Record

    Parameters
    ----------
    data
        Name of a .vcf file, or a list of `vcf.model._Record`
    chrom
        Chromosome name to be shown (only single sequence plot is supported).
    name
        Track name
    pad
        Size of paddings for each record. [`r.b - pad`..`r.b + pad`] is drawn.
    """
    if isinstance(data, str):
        data = load_vcf(data)
    if chrom is not None:
        data = list(filter(lambda x: x.CHROM == chrom, data))

    if len(data) == 0:
        return pl.lines(
            [(0, name, 0, name)], width=0, col=col, name=name, show_legend=False
        )

    return pl.lines(
        [(max(x.start - pad, 0), name, x.end + pad, name) for x in data],
        text=[
            (
                f"{x.CHROM} @ {x.POS}<br>{x.REF} -> {x.ALT}<br>{x.INFO}"
                if not single_sample
                else f"{x.CHROM} @ {x.POS}<br>{x.REF} -> {x.ALT}<br>{x.INFO}<br>{x.samples[0].sample}\t{x.samples[0].gt_type}\t{x.samples[0].gt_bases}"
            )
            for x in data
        ],
        width=width,
        col=col,
        name=name,
        show_legend=False,
        use_webgl=use_webgl,
    )
