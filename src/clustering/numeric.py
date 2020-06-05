from dataclasses import dataclass
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans, Birch
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import NMF
from logzero import logger
from .base import Clustering


@dataclass(repr=False, eq=False)
class ClusteringNumeric(Clustering):
    """Class for clustering of numerical data."""
    def calc_dist_mat(self, metric="euclidean"):
        """Required by hierarchical clustering and t-SNE plot."""
        self.c_dist_mat = pdist(self.data, metric=metric)
        self.s_dist_mat = squareform(self.c_dist_mat)

    ## ------------------------------------------ ##
    ##  Clustering with fixed number of clusters  ##
    ## ------------------------------------------ ##
    def kmeans(self, n_clusters):
        self.assignments = KMeans(n_clusters=n_clusters).fit_predict(self.data)

    def gmm(self, n_clusters):
        gmm = GaussianMixture(n_components=n_clusters, covariance_type="full")
        gmm.fit(self.data)
        self.assignments = gmm.predict(self.data)

    def birch(self, n_clusters):
        self.assignments = Birch(n_clusters=n_clusters).fit_predict(self.data)

    def nmf(self, n_clusters):   # TODO: debug
        nmf = NMF(n_components=n_clusters)
        W = nmf.fit_transform(self.data)
        H = nmf.components_
        print(W)
        print(H)

    ## ---------------------------------------------- ##
    ##  Clustering with arbitrary number of clusters  ##
    ## ---------------------------------------------- ##
    def _gmm_bic(self, n_clusters):
        gmm = GaussianMixture(n_components=n_clusters, covariance_type="full")
        gmm.fit(self.data)
        return gmm

    def gmm_bic(self, max_n_clusters=100):
        """Gaussian Mixture Model with model selection by BIC."""
        gmm_bics = [self._gmm_bic(i).bic(self.data) for i in range(1, max_n_clusters + 1)]
        logger.debug(gmm_bics)
        self.assignments = self._gmm_bic(np.argmin(gmm_bics) + 1).predict(self.data)
