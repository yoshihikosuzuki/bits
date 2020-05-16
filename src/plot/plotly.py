from typing import Any, Union, Optional, Sequence, Mapping, Tuple, List, Dict
from numbers import Number
from os.path import splitext
from collections import Counter
import numpy as np
import plotly.offline as py
import plotly.graph_objects as go
from plotly.basedatatypes import BaseTraceType


def make_line(x0: float, y0: float, x1: float, y1: float,
              width: float = 1, col: str = "black", layer: str = "above"):
    """Create a (non-interactive) line-shape object for Plotly.
    `layer` must be one of {"above" (default), "below"}.
    """
    assert layer in ("above", "below"), \
        "`layer` must be 'above' or 'below'"
    return dict(type="line", xref="x", yref="y", x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color=col, width=width), layer=layer)


def make_rect(x0: float, y0: float, x1: float, y1: float,
              xref: str = "x", yref: str = "y",
              fill_col: str = "grey", opacity: float = 1.,
              frame_width: float = 0, frame_col: Optional[str] = None,
              layer: str = "above"):
    """Create a (non-interactive) rectangle object for Plotly.
    `xref` must be one of {"x" (default), "paper"}. "paper" means `x0` and `x1` indicate
    horizontal relative positions of the entire plot (values are in [0, 1]).
    (Same goes for `yref`.)
    `layer` must be one of {"above" (default), "below"}.
    """
    assert xref in ("x", "paper") and yref in ("y", "paper"), \
        "`[x|y]ref` must be '[x|y]' or 'paper'"
    assert layer in ("above", "below"), \
        "`layer` must be 'above' or 'below'"
    return dict(type="rect", xref=xref, yref=yref, x0=x0, y0=y0, x1=x1, y1=y1,
                fillcolor=fill_col, opacity=opacity,
                line=dict(color=frame_col, width=frame_width), layer=layer)


def make_hist(data: Union[Sequence, Mapping[Any, int]],
              start: Optional[int] = None,
              end: Optional[int] = None,
              bin_size: Optional[int] = None,
              bin_num: int = 10,
              name: Optional[str] = None,
              show_legend: bool = False,
              use_histogram: bool = False):
    """Create a simple Trace object of a histogram for Plotly.
    Return not `go.Histogram` but `go.Bar` because the former is too inefficient.

    positional arguments:
      @ data : Raw data of numbers or counter of numbers

    optional arguments:
      @ start         : Start position of the plot range.
      @ end           : End position of the plot range.
      @ bin_size      : Size of each bin.
      @ bin_num       : Number of bins. Ignored if `bin_size` is set.
      @ name          : Display name of the trace in legend.
      @ show_legend   : Show this trace in legend.
      @ use_histogram : Force to use `go.Histogram`, not `go.Bar`.
    """
    if use_histogram:
        assert isinstance(data, Sequence), \
            "Only Sequence objects are supported if `use_histogram`."
        return go.Histogram(x=data, xbins=dict(start=start, end=end, size=bin_size),
                            name=name, showlegend=show_legend)

    if isinstance(data, Mapping) or not isinstance(next(iter(data)), Number):
        counter = Counter(data) if isinstance(data, Sequence) else data
        return go.Bar(**dict(zip(('x', 'y'), zip(*counter.items()))),
                      width=bin_size, name=name, showlegend=show_legend)

    if start is None:
        start = min(data)
    if end is None:
        end = max(data)
    if bin_size is not None:
        bin_num = -int(-(end - start) // bin_size)
    bin_size = (end - start) / bin_num

    counts, bin_edges = np.histogram(data, range=(start, end), bins=bin_num)

    return go.Bar(x=[(bin_edges[i] + bin_edges[i + 1]) / 2
                     for i in range(len(bin_edges) - 1)],
                  y=counts,
                  width=bin_size,
                  name=name, showlegend=show_legend)


def make_scatter(x: Sequence,
                 y: Sequence,
                 text: Optional[Sequence] = None,
                 text_pos: Optional[str] = None,
                 text_size: Optional[float] = None,
                 text_col: Optional[str] = None,
                 mode: str = "markers",
                 marker_size: float = 5,
                 line_width: float = 1,
                 col: Optional[str] = None,
                 col_scale: Optional[str] = None,
                 reverse_scale: bool = False,
                 show_scale: bool = True,
                 name: Optional[str] = None,
                 show_legend: bool = False):
    """Create a simple Trace object of a scatter plot for Plotly.

    positional arguments:
      @ x : Coordinates of data on x-axis.
      @ y : Coordinates of data on y-axis.

    optional arguments:
      @ text          : Texts for each data.
      @ text_pos      : Specify positions of `text`.
                        Format is: "[top|middle|bottom] [left|center|right]".
      @ text_size     : Size of `text`.
      @ text_col      : Color of `text`.
      @ mode          : Must be "markers", "lines", or "markers+lines".
      @ marker_size   : For "markers" mode.
      @ line_width    : For "lines" mode.
      @ col           : Color of markers and lines.
      @ col_scale     : Color scale for markers and lines.
      @ reverse_scale : Reverse `col_scale`.
      @ show_scale    : Show a scale bar for `col_scale`.
      @ name          : Display name of the trace in legend.
      @ show_legend   : Show this trace in legend.
    """
    assert len(x) == len(y), "`x` and `y` must have same size"
    if text is not None:
        assert len(x) == len(text), "`text` must have same size as data"
    assert mode in ("markers", "lines", "markers+lines"), \
        "`mode` must be one of {'markers', 'lines', 'markers+lines'}"
    return go.Scatter(x=x, y=y, text=text, mode=mode, name=name,
                      hoverinfo="text" if text is not None else None,
                      textposition=text_pos,
                      textfont=dict(size=text_size, color=text_col),
                      marker=(None if mode == "lines" else
                              dict(size=marker_size, color=col, colorscale=col_scale,
                                   reversescale=reverse_scale,
                                   showscale=show_scale if col_scale is not None else None)),
                      line=(None if mode == "markers" else
                            dict(width=line_width, color=col)),
                      showlegend=show_legend)


def make_layout(width: Optional[int] = None,
                height: Optional[int] = None,
                title: Optional[str] = None,
                x_title: Optional[str] = None,
                y_title: Optional[str] = None,
                x_range: Optional[Tuple[Optional[float], Optional[float]]] = None,
                y_range: Optional[Tuple[Optional[float], Optional[float]]] = None,
                x_grid: bool = True,
                y_grid: bool = True,
                x_zeroline: bool = True,
                y_zeroline: bool = True,
                x_reversed: bool = False,
                y_reversed: bool = False,
                shapes: Optional[List[Dict]] = None):
    """Create a simple Layout object for Plotly.

    optional arguments:
      @ width      : Width of the plot.
      @ height     : Height of the plot.
      @ title      : Title of the plot.
      @ x_title    : Title of x-axis of the plot.
      @ y_title    : Title of y-axis of the plot.
      @ x_range    : Range on x-axis to be drawn in the plot.
      @ y_range    : Range on y-axis to be drawn in the plot.
      @ x_grid     : Show grid on x-axis if True.
      @ y_grid     : Show grid on y-axis if True.
      @ x_zeroline : Show zeroline on x-axis if True.
      @ y_zeroline : Show zeroline on y-axis if True.
      @ x_reversed : Reverse x-axis if True.
      @ y_reversed : Reverse y-axis if True.
      @ shapes     : List of non-interactive shape objects.
    """
    return go.Layout(width=width, height=height, hovermode="closest", title=title,
                     xaxis=dict(title=x_title, showgrid=x_grid, zeroline=x_zeroline,
                                range=x_range,
                                autorange="reversed" if x_reversed else None),
                     yaxis=dict(title=y_title, showgrid=y_grid, zeroline=y_zeroline,
                                range=y_range,
                                autorange="reversed" if y_reversed else None),
                     shapes=shapes)


def show_plot(trace_list: Union[BaseTraceType, List[BaseTraceType]],
              layout: Optional[go.Layout] = None,
              out_fname: Optional[str] = None,
              include_plotlyjs: bool = False):
    """Plot a figure in Jupyter Notebook and/or to `out_fname`.

    positional arguments:
      @ trace_list : A trace or list of traces.

    optional arguments:
      @ layout           : A layout object.
      @ out_fname        : Output file name. Extension must be [png|jpeg|svg|html].
      @ include_plotlyjs : Embed JS codes for an independent plot if True.
    """
    out_type = None if out_fname is None else splitext(out_fname)[-1][1:]
    assert out_type is None or out_type in ("png", "jpeg", "svg", "html"), \
        "Unsupported output file type"

    fig = go.Figure(data=trace_list,
                    layout=make_layout() if layout is None else layout)
    py.iplot(fig, image=out_type if out_type != "html" else None,
             filename=None if out_fname is None else splitext(out_fname)[0])
    if out_type == "html":
        fig.write_html(file=out_fname, include_plotlyjs=include_plotlyjs)
