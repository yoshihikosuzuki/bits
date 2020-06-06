from copy import copy
from dataclasses import dataclass, InitVar
from typing import Sequence, List, Tuple
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
        """Compute a distance matrix among the input sequencs."""
        self.s_dist_mat = (calc_dist_mat_single(self.data, self.er) if n_core == 1
                           else calc_dist_mat_multi(self.data, self.er, n_core))
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
        self.cons_seqs = self._generate_consensus()
        logger.info(f"Initial consensus sequences:\n{self.cons_seqs}")
        self._merge_similar_clusters(th_merge)
        self._remove_small_clusters(th_noisy)
        self._sync_clusters(th_sync)
        logger.info(f"Final consensus sequences:\n{self.cons_seqs}")

    def _generate_consensus(self) -> List[ClusterCons]:
        """For each cluster, compute a consensus seuqnece among sequences
        belonging to the cluster."""
        # input sequences are ready to take consensus?
        is_aligned = not (self.cyclic or self.revcomp)
        return [ClusterCons(
            seq=consed.consensus(list(seqs) if is_aligned
                                 else [seq if i == 0
                                       else self.er.align(seqs[0], seq).b_aligned_seq
                                       for i, seq in enumerate(seqs)],
                                 seed_choice="median",
                                 n_iter=3),
            cluster_id=cluster_id,
            cluster_size=len(seqs))
            for cluster_id, seqs in self.clusters()]

    def _merge_similar_clusters(self, th_merge: float):
        """Merge clusters whose consensus sequences are similar until there
        exist no such clusters."""
        def __merge_similar_clusters() -> bool:
            flag_next = False
            for i, cons_i in enumerate(self.cons_seqs):
                for j, cons_j in enumerate(self.cons_seqs):
                    if i >= j:
                        continue
                    if self.er.align(cons_i.seq, cons_j.seq).diff < th_merge:
                        self.merge_cluster(cons_j.cluster_id,
                                           cons_i.cluster_id)
                        flag_next = True
            self.cons_seqs = self._generate_consensus()
            return flag_next

        while __merge_similar_clusters():
            pass

    def _remove_small_clusters(self, th_noisy: float):
        """Remove too small clusters."""
        min_size = max([1, self.N * th_noisy])
        self.cons_seqs = list(filter(lambda x: x.cluster_size <= min_size,
                                     self.cons_seqs))

    def _sync_clusters(self, th_sync: float):
        """Synchronize start positions of similar consensus sequences."""
        for i, cons_i in enumerate(self.cons_seqs):
            for j, cons_j in enumerate(self.cons_seqs):
                if i >= j:
                    continue
                aln = self.er.align(cons_i.seq, cons_j.seq)
                if aln.diff < th_sync:
                    cons_j.seq = aln.b_aligned_seq


def calc_dist_mat_single(seqs: Sequence[str],
                         er: EdlibRunner) -> np.ndarray:
    N = len(seqs)
    dist_mat = np.zeros((N, N), dtype="float32")
    for i, seq_i in seqs:
        for j, seq_j in seqs:
            if i >= j:
                continue
            dist_mat[i][j] = dist_mat[j][i] = er.align(seq_i, seq_j).diff
    return dist_mat


def calc_dist_mat_multi(seqs: Sequence[str],
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
