from dataclasses import dataclass, field
from typing import Any, List
from collections import Counter
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from BITS.plot.plotly import make_hist, make_scatter, make_layout, show_plot


@dataclass(repr=False, eq=False)
class Clustering:
    """Base class of clustering, offering basic and common variables and operations for clustering."""
    data        : Any   # array-like (List, pd.Series, np.ndarray((N, L))) data
    names       : List[str]  = None   # of each data; displayed in plots
    N           : int        = field(init=False)   # number of data; = len(data) or data.shape[0]
    L           : int        = field(init=False, default=None)   # number of features; = data.shape[1]
    assignments : np.ndarray = field(init=False)   # cluster assignment for each data; int type of length N
    s_dist_mat  : np.ndarray = field(init=False, default=None)   # square distance matrix
    c_dist_mat  : np.ndarray = field(init=False, default=None)   # condensed distance matrix
    cache       : dict       = field(default_factory=dict)   # store large intermediate data

    def __post_init__(self):
        if type(self.data) is list:
            self.data = pd.Series(self.data)
        self.N = self.data.shape[0]
        if len(self.data.shape) > 1:
            self.L = self.data.shape[1]
        if self.names is None:
            self.names = list(range(self.N))
        self.assignment = np.full(self.N, -1, dtype='int8')

    @property
    def n_clusters(self):
        return len(set(self.assignment))

    def cluster(self, cluster_id, return_where=False):
        """Return the data assigned to cluster <cluster_id>.
        If <return_where>, return the indices of the data instead of the data itself.
        """
        where = np.where(self.assignment == cluster_id)[0]
        return where if return_where else self.data[where]

    def clusters(self, return_where=False):
        """Generator of the clusters in order of cluster ID.
        If <return_where>, return the indices of the data instead of the data itself.
        """
        for cluster_id in sorted(list(set(self.assignment))):
            yield (cluster_id, self.cluster(cluster_id, return_where))

    def merge_cluster(self, cluster_from, cluster_to):
        """Change the assignments of the data in <cluster_from> as <cluster_to>."""
        self.assignment[self.cluster(cluster_from, return_where=True)] = cluster_to

    def hist_cluster_size(self, bin_size=10, width=18, height=10):
        """Histogram of the size of the clusters."""
        trace = make_hist(list(Counter(self.assignment).values()), bin_size=bin_size)
        layout = make_layout(width, height, x_title="Cluster size", y_title="Frequency")
        show_plot([trace], layout)

    def plot_dist_mat(self, variable_scale=False, show_scale=True, size=500, out_fname=None):
        """Draw a heatmap of <s_dist_mat>."""
        assert self.N <= 1000, "Too many data to plot"
        trace = go.Heatmap(z=self.s_dist_mat, colorscale="YlGnBu", showscale=show_scale,
                           zmin=None if variable_scale else 0., zmax=None if variable_scale else 1.)
        layout = make_layout(size, size, y_reversed=True)
        show_plot([trace], layout, out_fname=out_fname)

    def plot_tsne(self, size=700, marker_size=5, title=None, out_fname=None):
        """Draw t-SNE plot of the data. Requires <s_dist_mat>."""
        assert self.s_dist_mat is not None, "No distance matrix"
        coord = TSNE(n_components=2, metric='precomputed').fit_transform(self.s_dist_mat)
        trace = make_scatter(coord[:, 0].T, coord[:, 1].T,
                             text=[f"{self.names[i]}<br>cluster {self.assignment[i]}" for i in range(self.N)],
                             col=self.assignment, col_scale="Rainbow", show_scale=False,
                             marker_size=marker_size)
        layout = make_layout(size, size, title=title)
        show_plot([trace], layout, out_fname=out_fname)

    def _calc_dendrogram(self, method):
        if method not in self.cache:
            assert self.c_dist_mat is not None, "No condensed distance matrix"
            self.cache[method] = linkage(self.c_dist_mat, method=method)
            return self.cache[method]
        return self.cache[method]

    def show_dendrogram(self, method="ward", figsize=(18, 10)):
        plt.figure(figsize=figsize)
        dendrogram(self._calc_dendrogram(method))
        plt.show()

    def cluster_hierarchical(self, method="ward", criterion="distance", threshold=0.7, figsize=(18, 10)):
        """<method> = {single, complete, average, weighted, centroid, median, ward (default)}
        <criterion> = {inconsistent, distance (default), maxclust, monocrit, maxclust_monocrit}
        <threshold> = threshold in distance criterion
        """
        self.assignment = np.array(fcluster(self._calc_dendrogram(method), t=threshold, criterion=criterion))
