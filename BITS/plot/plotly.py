from os.path import splitext
import plotly.offline as py
import plotly.graph_objects as go


######################################
# --- Wrapper for shape objects. --- #
######################################

def make_line(x0, y0, x1, y1, col="black", width=1, layer="below"):
    """Create a (non-interactive) line-shape object for Plotly."""
    return dict(type="line", xref="x", yref="y", x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color=col, width=width), layer=layer)


def make_rect(x0, y0, x1, y1, xref="x", yref="y", fill_col="grey", opacity=1.,
              width=0, layer="below"):
    """Create a (non-interactive) rectangle object."""
    return dict(type="rect", xref=xref, yref=yref, x0=x0, y0=y0, x1=x1, y1=y1,
                fillcolor=fill_col, opacity=opacity, line=dict(width=width), layer=layer)


######################################
# --- Wrapper for Trace objects. --- #
######################################

def make_hist(x, start=None, end=None, bin_size=None, name=None, show_legend=True):
    """Create a Plotly trace object of Histogram plot."""
    return go.Histogram(x=x, xbins=dict(start=start, end=end, size=bin_size),
                        name=name, showlegend=show_legend)


def make_scatter(x, y, text=None, text_pos=None, text_size=None, text_col=None,
                 mode="markers", marker_size=5, col=None, col_scale=None, show_scale=True,
                 name=None, show_legend=True):
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

def make_layout(width=None, height=None, title=None, x_title=None, y_title=None,
                x_range=None, y_range=None, x_grid=True, y_grid=True,
                x_zeroline=True, y_zeroline=True, x_reversed=False, y_reversed=False,
                shapes=None):
    """Create a Plotly Layout object."""
    return go.Layout(width=width, height=height, hovermode="closest", title=title,
                     xaxis=dict(title=x_title, showgrid=x_grid, zeroline=x_zeroline, range=x_range,
                                autorange="reversed" if x_reversed else None),
                     yaxis=dict(title=y_title, showgrid=y_grid, zeroline=y_zeroline, range=y_range,
                                autorange="reversed" if y_reversed else None),
                     shapes=shapes)


#################################
# --- Wrapper for Plotting. --- #
#################################

def show_plot(trace_list, layout=None, out_fname=None, include_plotlyjs=False):
    """Plot a figure in a Notebook or to <out_fname>.
    File extension of <out_fname> can be "png", "jpeg", "svg" or "html".
    <include_plotlyjs> = True, False or "cdn". True increases the file size by 3MB.
    """
    out_type = None if out_fname is None else splitext(out_fname)[-1][1:]
    assert out_type is None or out_type in ("png", "jpeg", "svg", "html"), "Not supported file type"

    fig = go.Figure(data=trace_list, layout=make_layout() if layout is None else layout)
    py.iplot(fig, image=out_type if out_type != "html" else None,
             filename=None if out_fname is None else splitext(out_fname)[0])
    if out_type == "html":
        fig.write_html(file=out_fname, include_plotlyjs=include_plotlyjs)
