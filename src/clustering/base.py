from dataclasses import dataclass, field
from typing import Optional, Sequence, List, Tuple
from collections import Counter
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from logzero import logger
from BITS.plot.plotly import make_hist, make_scatter, make_layout, show_plot


@dataclass(repr=False, eq=False)
class Clustering:
    """Base class of clustering, offering common variables and operations.

    positional variables:
      @ data  : Array of data to be clustered.

    optional variables:
      @ names : Display name for each data.

    uninitialized variables:
      @ N           : Number of data.
      @ assignments : Cluster assignment for each data.
      @ s_dist_mat  : Square distance matrix among data.
      @ c_dist_mat  : Condensed distance matrix among data.
      @ cache       : Internal cache.
    """
    data: Sequence
    names: Optional[Sequence[str]] = None
    N: int = field(init=False)
    assignments: np.ndarray = field(init=False)
    s_dist_mat: np.ndarray = field(init=False, default=None)
    c_dist_mat: np.ndarray = field(init=False, default=None)
    cache: dict = field(default_factory=dict)

    def __post_init__(self):
        self.N = len(self.data)
        logger.info(f"Input data size = {self.N}")
        if isinstance(self.data, list):
            self.data = np.array(self.data)
        if self.names is None:
            self.names = list(map(lambda i: f"data[{i}]", range(self.N)))
        self.assignments = np.full(self.N, -1, dtype=np.int16)

    @property
    def n_clusters(self) -> int:
        return len(set(self.assignments))

    def cluster(self,
                cluster_id: int,
                return_where: bool = False) -> np.ndarray:
        """Return the data assigned to cluster `cluster_id`.
        If `return_where`, return the indices of the data instead."""
        where = np.where(self.assignments == cluster_id)[0]
        return where if return_where else self.data[where]

    def clusters(self,
                 return_where: bool = False) -> List[Tuple[int, np.ndarray]]:
        """Generator of the clusters in order of cluster ID.
        If `return_where`, return the indices of the data instead."""
        return [(cluster_id, self.cluster(cluster_id, return_where))
                for cluster_id in sorted(set(self.assignments))]

    def merge_cluster(self,
                      cluster_from: int,
                      cluster_to: int):
        """Change the assignments of the data in `cluster_from` as `cluster_to`."""
        indices = self.cluster(cluster_from, return_where=True)
        self.assignments[indices] = cluster_to

    def hist_cluster_size(self,
                          bin_size: int = 1,
                          size: int = 500):
        """Histogram of the size of the clusters."""
        show_plot(make_hist(list(Counter(self.assignments).values()),
                            bin_size=bin_size),
                  make_layout(size,
                              size,
                              x_title="Cluster size",
                              y_title="Frequency"))

    def plot_dist_mat(self,
                      variable_scale: bool = False,
                      show_scale: bool = True,
                      size: int = 500):
        """Draw a heatmap of `s_dist_mat`."""
        assert self.s_dist_mat is not None, \
            "Square distance matrix is required"
        if self.N > 1000:
            logger.warn(f"Data size {self.N} is too large for heatmap."
                        f"Continue? [y/N]")
            if input() != "y":
                return
        show_plot(go.Heatmap(z=self.s_dist_mat,
                             colorscale="Blues",
                             reversescale=True,
                             showscale=show_scale,
                             zmin=None if variable_scale else 0.,
                             zmax=None if variable_scale else 1.),
                  make_layout(size,
                              size,
                              y_reversed=True,
                              margin=dict(l=10, r=10, t=50, b=10)))

    def plot_tsne(self,
                  size: int = 700,
                  marker_size: int = 5):
        """Draw t-SNE plot of the data."""
        assert self.s_dist_mat is not None, \
            "Square distance matrix is required"
        coord = (TSNE(n_components=2, metric="precomputed")
                 .fit_transform(self.s_dist_mat))
        show_plot(make_scatter(coord[:, 0].T,
                               coord[:, 1].T,
                               text=[f"{self.names[i]}<br>"
                                     f"cluster {self.assignments[i]}"
                                     for i in range(self.N)],
                               col=self.assignments,
                               col_scale="Rainbow",
                               show_scale=False,
                               marker_size=marker_size),
                  make_layout(width=size,
                              height=size,
                              x_grid=False,
                              y_grid=False,
                              x_zeroline=False,
                              y_zeroline=False,
                              x_show_tick_labels=False,
                              y_show_tick_labels=False,
                              margin=dict(l=10, r=10, t=50, b=10)))

    def cluster_hierarchical(self,
                             threshold: float,
                             method: str = "ward",
                             criterion: str = "distance"):
        """Run a hierarchical clustering.

        positional arguments:
          @ threshold : Threshold for `criterion`.

        optional arguments:
          @ method    : Must be one of  {"single", "complete", "average",
                        "weighted", "centroid", "median", "ward"}.
          @ criterion : Must be one of {"inconsistent", "distance", "maxclust",
                        "monocrit", "maxclust_monocrit"}.
        """
        self.assignments = np.array(fcluster(self._calc_dendrogram(method),
                                             t=threshold,
                                             criterion=criterion))
        logger.info(f"{self.n_clusters} clusters generated")

    def _calc_dendrogram(self,
                         method: str) -> np.ndarray:
        if method not in self.cache:
            assert self.c_dist_mat is not None, \
                "Condensed distance matrix is required"
            self.cache[method] = linkage(self.c_dist_mat, method=method)
        return self.cache[method]

    def show_dendrogram(self,
                        method: str = "ward",
                        figsize: Tuple[int, int] = (18, 10)):
        """Draw a dendrogram of the hierarchical clustering."""
        plt.figure(figsize=figsize)
        dendrogram(self._calc_dendrogram(method))
        plt.show()
