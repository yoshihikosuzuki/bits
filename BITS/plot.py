import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
from logzero import logger


def make_line(x0, y0, x1, y1, col='black', width=1, layer='below'):
    """
    For Plotly.
    Create a line-shape object for Plotly.
    """

    return {'type': 'line', 'xref': 'x', 'yref': 'y',
            'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
            'line': {'color': col, 'width': width},
            'layer': layer}


def plot_scatter(x, y, text=None, marker_size=5, width=1000, height=1000, x_title=None, y_title=None):
    trace = go.Scatter(x=x,
                       y=y,
                       mode="markers",
                       marker=dict(size=marker_size))
    if text is not None:
        trace["text"] = text
        trace["hoverinfo"] = "text"

    layout = go.Layout(width=width,
                       height=height,
                       hovermode='closest')
    if x_title is not None:
        layout["xaxis"] = dict(title=x_title)
    if y_title is not None:
        layout["yaxis"] = dict(title=y_title)

    py.iplot(go.Figure(data=[trace], layout=layout))
