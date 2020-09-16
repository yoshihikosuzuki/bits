from typing import Any, Union, Optional, Sequence, Mapping, Tuple, List, Dict
from numbers import Number
from collections import Counter
import numpy as np
import plotly.graph_objects as go
from plotly.basedatatypes import BaseTraceType


def make_line(x0: float, y0: float, x1: float, y1: float,
              width: float = 1,
              col: str = "black",
              layer: str = "above") -> Dict:
    """Create a (non-interactive) line-shape object for Plotly.
    `layer` must be one of {"above" (default), "below"}.
    """
    assert layer in ("above", "below"), \
        "`layer` must be 'above' or 'below'"
    return dict(type="line",
                xref="x",
                yref="y",
                x0=x0,
                y0=y0,
                x1=x1,
                y1=y1,
                line=dict(color=col,
                          width=width),
                layer=layer)


def make_lines(x0y0x1y1s: List[Tuple[int, int, int, int]],
               width: float = 1,
               col: str = "black",
               name: Optional[str] = None,
               show_legend: bool = False,
               use_webgl: bool = False) -> go.Scatter:
    """Scatter trace object is lighter than shape object. However, only single
    width and color is allowed."""
    return make_scatter(x=[x for x0, _, x1, _ in x0y0x1y1s
                           for x in (x0, x1, None)],
                        y=[y for _, y0, _, y1 in x0y0x1y1s
                           for y in (y0, y1, None)],
                        mode="lines",
                        line_width=width,
                        col=col,
                        name=name,
                        show_legend=show_legend,
                        use_webgl=use_webgl)


def make_rect(x0: float, y0: float, x1: float, y1: float,
              xref: str = "x",
              yref: str = "y",
              fill_col: str = "grey",
              opacity: float = 1.,
              frame_width: float = 0,
              frame_col: Optional[str] = None,
              layer: str = "above") -> Dict:
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
    return dict(type="rect",
                xref=xref,
                yref=yref,
                x0=x0,
                y0=y0,
                x1=x1,
                y1=y1,
                fillcolor=fill_col,
                opacity=opacity,
                line=dict(color=frame_col,
                          width=frame_width),
                layer=layer)


def make_hist(data: Union[Sequence, Mapping[Any, int]],
              start: Optional[int] = None,
              end: Optional[int] = None,
              bin_size: Optional[int] = None,
              bin_num: int = 10,
              name: Optional[str] = None,
              show_legend: bool = False,
              use_histogram: bool = False) -> Union[go.Bar, go.Histogram]:
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
        return go.Histogram(x=data,
                            xbins=dict(start=start,
                                       end=end,
                                       size=bin_size),
                            name=name,
                            showlegend=show_legend)

    if isinstance(data, Mapping) or not isinstance(next(iter(data)), Number):
        counter = Counter(data) if isinstance(data, Sequence) else data
        return go.Bar(**dict(zip(('x', 'y'),
                                 zip(*counter.items()))),
                      width=bin_size,
                      name=name,
                      showlegend=show_legend)

    if len(set(data)) == 1:
        return go.Bar(x=[data[0]],
                      y=[len(data)],
                      width=bin_size,
                      name=name,
                      showlegend=show_legend)

    if start is None:
        start = min(data)
    if end is None:
        end = max(data)
    if bin_size is not None:
        bin_num = -int(-(end - start) // bin_size)
    bin_size = (end - start) / bin_num

    counts, bin_edges = np.histogram(data,
                                     range=(start, end),
                                     bins=bin_num)

    return go.Bar(x=[(bin_edges[i] + bin_edges[i + 1]) / 2
                     for i in range(len(bin_edges) - 1)],
                  y=counts,
                  width=bin_size,
                  name=name,
                  showlegend=show_legend)


def make_scatter(x: Sequence,
                 y: Sequence,
                 text: Optional[Sequence] = None,
                 text_pos: Optional[str] = None,
                 text_size: Optional[float] = None,
                 text_col: Optional[str] = None,
                 mode: str = "markers",
                 marker_size: float = 5,
                 marker_width: int = None,
                 line_width: float = 1,
                 col: Optional[Union[str, Sequence[str]]] = None,
                 col_range: Optional[Tuple[float, float]] = None,
                 col_scale: Optional[str] = None,
                 reverse_scale: bool = False,
                 show_scale: bool = True,
                 name: Optional[str] = None,
                 show_legend: bool = False,
                 use_webgl: bool = False) -> go.Scatter:
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
      @ mode          : "markers", "lines", "markers+lines", "text", etc.
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
    return (go.Scattergl if use_webgl else go.Scatter)(
        x=x,
        y=y,
        text=text,
        mode=mode,
        name=name,
        hoverinfo="text" if text is not None else None,
        textposition=text_pos,
        textfont=dict(size=text_size,
                      color=text_col),
        marker=(None if mode == "lines" else
                dict(size=marker_size,
                     color=col,
                     cmin=col_range[0] if col_range is not None else None,
                     cmax=col_range[1] if col_range is not None else None,
                     colorscale=col_scale,
                     reversescale=reverse_scale,
                     showscale=show_scale if col_scale is not None else None,
                     line=dict(width=marker_width))),
        line=(None if mode == "markers" else
              dict(width=line_width,
                   color=col)),
        showlegend=show_legend)


def make_layout(width: Optional[int] = None,
                height: Optional[int] = None,
                font: str = "Arial",
                font_size_title: int = 20,
                font_size_axis_title: int = 18,
                font_size_axis_tick: int = 15,
                font_size_legend: int = 15,
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
                x_show_tick_labels: bool = True,
                y_show_tick_labels: bool = True,
                anchor_axes: bool = False,
                margin: Dict = dict(l=80, r=80, t=100, b=80),
                shapes: Optional[List[Dict]] = None) -> go.Layout:
    """Create a simple Layout object for Plotly.

    optional arguments:
      @ width       : Width of the plot.
      @ height      : Height of the plot.
      @ font        : Font of (axis) titles.
      @ font_size_[title|axis_title|axis_tick]
                    : Font size of title/[x|y]_title/[x|y]_tick.
      @ title       : Title of the plot.
      @ x_title     : Title of x-axis of the plot.
      @ y_title     : Title of y-axis of the plot.
      @ x_range     : Range on x-axis to be drawn in the plot.
      @ y_range     : Range on y-axis to be drawn in the plot.
      @ x_grid      : Show grid on x-axis if True.
      @ y_grid      : Show grid on y-axis if True.
      @ x_zeroline  : Show zeroline on x-axis if True.
      @ y_zeroline  : Show zeroline on y-axis if True.
      @ x_reversed  : Reverse x-axis if True.
      @ y_reversed  : Reverse y-axis if True.
      @ anchor_axes : Use same scale for both x/y axes.
      @ margin      : Margin of the plot.
      @ shapes      : List of non-interactive shape objects.
    """
    return go.Layout(
        width=width,
        height=height,
        hovermode="closest",
        title=dict(text=title,
                   font=dict(family=font,
                             size=font_size_title,
                             color="black")),
        xaxis=dict(title=dict(text=x_title,
                              font=dict(family=font,
                                        size=font_size_axis_title,
                                        color="black")),
                   tickfont=dict(family="Arial",
                                 size=font_size_axis_tick,
                                 color="black"),
                   showgrid=x_grid,
                   zeroline=x_zeroline,
                   showticklabels=x_show_tick_labels,
                   range=x_range,
                   autorange="reversed" if x_reversed else None),
        yaxis=dict(title=dict(text=y_title,
                              font=dict(family=font,
                                        size=font_size_axis_title,
                                        color="black")),
                   tickfont=dict(family="Arial",
                                 size=font_size_axis_tick,
                                 color="black"),
                   showgrid=y_grid,
                   zeroline=y_zeroline,
                   showticklabels=y_show_tick_labels,
                   range=y_range,
                   autorange="reversed" if y_reversed else None,
                   scaleanchor="x" if anchor_axes else None),
        legend=dict(font=dict(family=font,
                              size=font_size_legend,
                              color="black")),
        margin=margin,
        shapes=shapes)


def show_plot(traces: Union[BaseTraceType, List[BaseTraceType]],
              layout: Optional[go.Layout] = None,
              download_as: str = "svg",
              out_html: Optional[str] = None,
              embed_plotlyjs: bool = True):
    """Plot a figure in Jupyter Notebook.

    positional arguments:
      @ traces : A trace or list of traces.

    optional arguments:
      @ layout         : A layout object.
      @ download_as    : File format of the "Download plot" buttion in the plot.
      @ out_html       : Output the plot to an html file.
      @ embed_plotlyjs : If True, embed plotly.js codes (~3 MB) in `out_html`.
    """
    assert download_as in ("png", "jpeg", "svg"), \
        f"Unsupported output file type: {download_as}"
    fig = go.Figure(data=traces,
                    layout=layout if layout is not None else make_layout())
    fig.show(config=dict(toImageButtonOptions=dict(format=download_as)))
    if out_html is not None:
        fig.write_html(file=out_html, include_plotlyjs=embed_plotlyjs)
