from dataclasses import dataclass
import numpy as np
import pandas as pd
from .numeric import ClusteringNumeric


@dataclass(repr=False, eq=False)
class ClusteringVarMat(ClusteringNumeric):   # TODO: obsolete class?
    def calc_dist_mat(self):
        super().calc_dist_mat(metric="hamming")

    def generate_consensus(self, how="representative"):
        """
        Compute consensus for each cluster.

        <how> == "representative": majority vote,
              == "centroid": average (mainly for debug)
        """

        assert how in ("representative", "centroid"), "<how> must be 'representative' or 'centroid'"

        ret = {}
        index = 0
        for cluster_id, cluster_data in super().clusters():
            n_cluster = cluster_data.shape[0]
            ret[index] = (cluster_id,
                          n_cluster,
                          np.where(np.sum(cluster_data, axis=0) >= n_cluster // 2, 1, 0) if how == "representative"
                          else np.sum(cluster_data, axis=0) / n_cluster)

        return pd.DataFrame.from_dict(ret, orient="index",
                                      columns=("cluster_id", "cluster_size", "representative"))
