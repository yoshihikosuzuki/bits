from copy import copy
from dataclasses import dataclass, InitVar
from typing import List, Tuple
import numpy as np
from scipy.spatial.distance import squareform
from logzero import logger
import consed
from BITS.seq.io import SeqRecord
from BITS.seq.align import EdlibRunner
from BITS.util.proc import NoDaemonPool
from .base import Clustering


@dataclass
class ClusterCons(SeqRecord):
    cluster_id: int
    cluster_size: int

    def __repr__(self):
        return self._order_repr(["cluster_id", "cluster_size", "length", "seq"])


@dataclass(repr=False, eq=False)
class ClusteringSeq(Clustering):
    """Class for clustering of DNA sequences."""
    revcomp: InitVar[bool] = True
    cyclic: InitVar[bool] = False

    def __post_init__(self,
                      revcomp: bool,
                      cyclic: bool):
        super().__post_init__()
        self.er = EdlibRunner("global", revcomp=revcomp, cyclic=cyclic)

    def calc_dist_mat(self,
                      n_core: int = 1):
        self.s_dist_mat = calc_dist_mat(self.data, self.er, n_core)
        self.c_dist_mat = squareform(self.s_dist_mat)

    def generate_consensus(self,
                           th_merge: float = 0.05,
                           th_noisy: float = 0.01,
                           th_sync: float = 0.3):
        """
        1. Compute consensus sequences using Consed for each cluster.
        2. Merge every two consensus sequences with dissimilarity < `th_merge`.
        3. Remove every cluster whose size is < `N * th_noisy`.
        4. Synchronize consensus sequences with dissimilarity < `th_sync`.
        """
        assert np.max(self.assignments) >= 0, "Do clustering prior to this"
        # Initial consensus sequences
        self.cons_seqs = self._generate_consensus()
        logger.info(f"Initial consensus sequences:\n{self.cons_seqs}")
        # Merge too similar clusters
        self._merge_clusters(th_merge)
        # Remove too small clusters
        self.cons_seqs = list(filter(lambda x:
                                     x.cluster_size < max([2, self.N * th_noisy])))
        # Synchronize the consensus units
        n_cons = len(self.cons_seqs)
        for i in range(n_cons - 1):
            for j in range(i + 1, n_cons):
                aln = self.er.align(self.cons_seqs[i].seq,
                                    self.cons_seqs[j].seq)
                if aln.diff < th_sync:
                    self.cons_seqs[j].seq = aln.b_aligned_seq
        logger.info(f"Final consensus sequences:\n{self.cons_seqs}")

    def _generate_consensus(self) -> List[ClusterCons]:
        """For each cluster, compute a consensus seuqnece among sequences
        belonging to the cluster."""
        synchronized = not (self.cyclic or self.revcomp)
        return [ClusterCons(
            seq=consed.consensus(list(seqs) if synchronized
                                 else [seq if i == 0
                                       else self.er.align(seqs[0], seq).b_aligned_seq
                                       for i, seq in enumerate(seqs)],
                                 seed_choice="median",
                                 n_iter=3),
            cluster_id=cluster_id,
            cluster_size=seqs.shape[0])
            for cluster_id, seqs in self.clusters()]

    def _merge_clusters(self,
                        th_merge: float):
        """Merge clusters whose consensus sequences are similar until there
        exist no such clusters."""
        flag_next = True
        while flag_next:
            n_cons = len(self.cons_seqs)
            for i in range(n_cons - 1):
                for j in range(i + 1, n_cons):
                    if self.er.align(self.cons_seqs[i].seq,
                                     self.cons_seqs[j].seq).diff < th_merge:
                        self.merge_cluster(self.cons_seqs[j].cluster_id,
                                           self.cons_seqs[i].cluster_id)
                        flag_next = True
            self.cons_seqs = self._generate_consensus()


def calc_dist_mat(seqs: List[str],
                  er: EdlibRunner,
                  n_core: int = 1) -> np.ndarray:
    global global_seqs
    global_seqs = copy(seqs)
    N = len(seqs)
    rows = [i // 2 if i % 2 == 0 else N - 1 - i // 2 for i in range(N)]
    n_unit = -(-N // n_core)
    dist_mat = np.zeros((N, N), dtype="float32")
    with NoDaemonPool(n_core) as pool:
        for ret in pool.starmap(calc_dist_arrays,
                                [(rows[i * n_unit:(i + 1) * n_unit], er)
                                 for i in range(n_core)]):
            for row, dist_array in ret:
                dist_mat[row, row + 1:] = dist_mat[row + 1:, row] = dist_array
    return dist_mat


def calc_dist_arrays(rows: List[int],
                     er: EdlibRunner) -> List[Tuple[int, np.ndarray]]:
    return [(row,
             np.array([er.align(global_seqs[row], global_seqs[col]).diff
                       for col in range(row + 1, len(global_seqs))],
                      dtype="float32"))
            for row in rows]
