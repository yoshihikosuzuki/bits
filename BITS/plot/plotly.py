import plotly.offline as py
import plotly.graph_objs as go


def make_line(x0, y0, x1, y1, col="black", width=1, layer="below"):
    """Create a line-shape object for Plotly."""
    return dict(type="line", xref="x", yref="y", x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color=col, width=width), layer=layer)


######################################
# --- Wrappers for Trace object. --- #
######################################

def make_hist(x, start=None, end=None, bin_size=None):
    """Create a Plotly trace object of Histogram plot."""
    return go.Histogram(x=x, xbins=dict(start=start, end=end, size=bin_size))


def make_scatter(x, y, text=None, text_pos=None, text_size=None, text_col=None, mode="markers",
                 col=None, col_scale=None, show_scale=True, marker_size=5, name=None, show_legend=True):
    """Create a Plotly trace object of Scatter plot."""
    return go.Scatter(x=x, y=y, text=text, mode=mode, name=name,
                      hoverinfo="text" if text is not None else None,
                      textposition=text_pos, textfont=dict(size=text_size, color=text_col),
                      marker=dict(size=marker_size, color=col, colorscale=col_scale,
                                  showscale=show_scale if col_scale is not None else None),
                      showlegend=show_legend)


######################################
# --- Wrapper for Layout object. --- #
######################################

def make_layout(width, height, title=None, x_title=None, y_title=None,
                x_range=None, y_range=None, x_grid=True, y_grid=True,
                x_reversed=False, y_reversed=False, shapes=None):
    """Create a Plotly Layout object."""
    return go.Layout(width=width, height=height, hovermode="closest", title=title,
                     xaxis=dict(title=x_title, showgrid=x_grid, range=x_range,
                                autorange="reversed" if x_reversed else None),
                     yaxis=dict(title=y_title, showgrid=y_grid, range=y_range,
                                autorange="reversed" if y_reversed else None),
                     shapes=shapes)


#################################
# --- Wrapper for Plotting. --- #
#################################

def show_plot(trace_list, layout, out_fname=None):
    """Plot a figure in a Notebook or to <out_fname>."""
    fig = go.Figure(data=trace_list, layout=layout)
    py.iplot(fig) if out_fname is None else py.iplot(fig, filename=out_fname)


