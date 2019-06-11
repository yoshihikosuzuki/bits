import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
from logzero import logger


def generate_layout(width,
                    height,
                    title=None,
                    x_title=None,
                    y_title=None):
    """Utility for creating the simplest Plotly Layout object."""
    layout = go.Layout(width=width,
                       height=height,
                       hovermode='closest')
    if title is not None:
        layout["title"] = title
    if x_title is not None:
        layout["xaxis"] = dict(title=x_title)
    if y_title is not None:
        layout["yaxis"] = dict(title=y_title)
    return layout


def show_plot(trace_list, layout, out_fname=None):
    """Utility for plotting a figure regardless of <out_fname>."""
    fig = go.Figure(data=trace_list, layout=layout)
    py.iplot(fig) if out_fname is None else py.iplot(fig, filename=out_fname)


def generate_scatter(x, y, text=None, marker_size=5, name=None):
    """Utility for creating the simplest Plotly Scatter object."""
    trace = go.Scatter(x=x,
                       y=y,
                       mode="markers",
                       marker=dict(size=marker_size))
    if text is not None:
        trace["text"] = text
        trace["hoverinfo"] = "text"
    if name is not None:
        trace["name"] = name
    return trace


def plot_scatter(x,
                 y,
                 text=None,
                 marker_size=5,
                 width=700,
                 height=700,
                 title=None,
                 x_title=None,
                 y_title=None,
                 out_fname=None):
    """Utility for plotting the simplest interactive scatter plot."""
    trace = generate_scatter(x, y, text=text, marker_size=marker_size)
    layout = generate_layout(width, height, title=title, x_title=x_title, y_title=y_title)
    show_plot([trace], layout, out_fname=out_fname)


def make_line(x0, y0, x1, y1, col='black', width=1, layer='below'):
    """Create a line-shape object for Plotly."""
    return {'type': 'line', 'xref': 'x', 'yref': 'y',
            'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
            'line': {'color': col, 'width': width},
            'layer': layer}
