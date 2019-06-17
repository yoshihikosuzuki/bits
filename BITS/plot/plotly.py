import plotly.offline as py
import plotly.graph_objs as go


def make_layout(width, height, title=None, x_title=None, y_title=None,
                x_reversed=False, y_reversed=False):
    """Create a Plotly Layout object."""
    return go.Layout(width=width, height=height, hovermode="closest", title=title,
                     xaxis=dict(title=x_title, autorange="reversed" if x_reversed else None),
                     yaxis=dict(title=y_title, autorange="reversed" if y_reversed else None))


def show_plot(trace_list, layout, out_fname=None):
    """Plot a figure in a Notebook or to <out_fname>."""
    fig = go.Figure(data=trace_list, layout=layout)
    py.iplot(fig) if out_fname is None else py.iplot(fig, filename=out_fname)


def make_line(x0, y0, x1, y1, col='black', width=1, layer='below'):
    """Create a line-shape object for Plotly."""
    return dict(type="line", xref="x", yref="y", x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color=col, width=width), layer=layer)


def make_scatter(x, y, text=None, col=None, col_scale=None, show_scale=True, marker_size=5, name=None):
    """Create a Plotly trace object of Scatter plot."""
    return go.Scatter(x=x, y=y, text=text, mode="markers", name=name,
                      hoverinfo="text" if text is not None else None,
                      marker=dict(size=marker_size, color=col, colorscale=col_scale,
                                  showscale=show_scale if col_scale is not None else None))


def plot_scatter(x, y, text=None, marker_size=5, width=700, height=700,
                 title=None, x_title=None, y_title=None, out_fname=None):
    """Plot the simplest interactive scatter plot."""
    trace = make_scatter(x, y, text=text, marker_size=marker_size)
    layout = make_layout(width, height, title=title, x_title=x_title, y_title=y_title)
    show_plot([trace], layout, out_fname=out_fname)


