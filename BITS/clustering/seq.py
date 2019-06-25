import argparse
from os.path import join
from dataclasses import dataclass, InitVar
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from logzero import logger
import consed
from BITS.seq.align import EdlibRunner
from BITS.util.proc import run_command, NoDaemonPool
from BITS.util.log import print_log
from BITS.util.io import save_pickle, load_pickle
from BITS.util.scheduler import Scheduler
from .base import Clustering


def __calc_dist_array(i, data, runner):
    '''Compute row i vs columns (i+1) to N in the distance matrix.'''
    logger.debug(f'Computing row {i}')
    return (i, np.array([runner.align(data[i], data[j]).diff for j in range(i + 1, data.shape[0])],
                        dtype='float32'))


def _calc_dist_array(rows, data, runner):
    return [__calc_dist_array(row, data, runner) for row in rows]


@dataclass(repr=False, eq=False)
class ClusteringSeq(Clustering):
    '''Class for clustering of DNA sequences.'''
    revcomp : InitVar[bool] = True   # allow reverse complement of the target sequence if True
    cyclic  : InitVar[bool] = False   # allow cyclic alignment between two sequences if True

    def __post_init__(self, revcomp, cyclic):
        super().__post_init__()
        self.runner = EdlibRunner('global', revcomp=revcomp, cyclic=cyclic)

    def _calc_dist_mat(self, row_indices, n_core):
        '''Compute distance @ row vs [row + 1, ..., self.N] for each row in <row_indices> in parallel.'''
        unit_n = -(-len(row_indices) // n_core)
        dist_arrays = []
        with NoDaemonPool(n_core) as pool:
            for dist_array in pool.starmap(_calc_dist_array,
                                           [(row_indices[i * unit_n:(i + 1) * unit_n],
                                             self.data,
                                             self.runner)
                                            for i in range(n_core)]):
                dist_arrays += dist_array
        return dist_arrays

    @print_log('distance matrix calculation')
    def calc_dist_mat(self, n_core=1, n_distribute=1, dir_name='.', out_prefix='clustering',
                      job_scheduler='sge', submit_command='qsub'):
        '''Compute all-vs-all distance matrix in parallel.'''
        if n_core * n_distribute > self.N:   # no need for parallelization
            n_core = n_distribute = 1

        # Reorder the rows in the distance matrix so that average cumulative size becomes uniform
        rows = [int(i / 2) if i % 2 == 0  else self.N - 2 - int((i - 1) / 2)
                for i in range(self.N - 1)]   # data[self.N] will not be query, so excluded

        if n_distribute == 1:
            dist_arrays = self._calc_dist_mat(rows, n_core)
        else:
            # Output this clustering object so that each job can access the data
            obj_fname = join(dir_name, f'{out_prefix}_obj.pkl')
            save_pickle(self, obj_fname)

            # Split into and submit jobs
            s = Scheduler(job_scheduler, submit_command)
            jids = []
            unit_n = -(-len(rows) // n_distribute)
            for i in range(n_distribute):
                index = str(i + 1).zfill(int(np.log10(n_distribute) + 1))
                rows_part = rows[i * unit_n:(i + 1) * unit_n]
                rows_fname = join(dir_name, f'{out_prefix}_rows.{index}.pkl')
                out_fname = join(dir_name, f'{out_prefix}.{index}.pkl')
                script_fname = join(dir_name, f'{out_prefix}.sh.{index}')
                save_pickle(rows_part, rows_fname)
                script = (f'python -m BITS.clustering.seq {obj_fname} {rows_fname} {out_fname} {n_core}')

                jids.append(s.submit(script,
                                     script_fname,
                                     job_name='calc_dist_mat',
                                     n_core=n_core))

            # Merge the results
            s.submit('sleep 1s',
                     join(dir_name, 'gather.sh'),
                     job_name='gather_dist_mat',
                     depend=jids,
                     wait=True)

            dist_arrays = []
            for fname in run_command(f'find {dir_name} -maxdepth 1 -name "{out_prefix}.*.pkl"').strip().split('\n'):
                dist_arrays += load_pickle(fname)

        # Store the arrays into matrix
        self.s_dist_mat = np.zeros((self.N, self.N), dtype='float32')
        for i, dist_array in dist_arrays:
            self.s_dist_mat[i, i + 1:] = self.s_dist_mat[i + 1:, i] = dist_array
        self.c_dist_mat = squareform(self.s_dist_mat)

    def _generate_consensus(self):
        synchronized = not (self.cyclic or self.revcomp)
        ret = {}
        index = 0
        for cluster_id, seqs in self.clusters():
            if len(seqs) <= 1:
                logger.info(f'Skip isolated cluster {cluster_id}.')
                continue
            # TODO: better way to choose initial seed sequence instead of the first one?
            cons_seq = consed.consensus(list(seqs) if synchronized
                                        else [seq if i == 0
                                              else self.runner.align(seqs.iloc[0], seq).mapped_seq(seq)
                                              for i, seq in enumerate(seqs)],
                                        n_iter=3)
            if cons_seq == '':
                logger.warn(f'Consed failed: cluster {cluster_id} ({len(seqs)} seqs)')
                continue
            ret[index] = (cluster_id, seqs.shape[0], len(cons_seq), cons_seq)
            index += 1
        return pd.DataFrame.from_dict(ret, orient='index',
                                      columns=('cluster_id', 'cluster_size', 'length', 'sequence'))

    def _merge_clusters(self, th_merge):
        flag_next = False
        n_cons = self.cons_seqs.shape[0]   # != <self.n_clusters> due to Consed error
        for i in range(n_cons - 1):
            for j in range(i + 1, n_cons):
                if self.runner.align(self.cons_seqs['sequence'].iloc[i],
                                     self.cons_seqs['sequence'].iloc[j]).diff < th_merge:
                    self.merge_cluster(self.cons_seqs['cluster_id'].iloc[j],
                                       self.cons_seqs['cluster_id'].iloc[i])
                    flag_next = True
        self.cons_seqs = self._generate_consensus()
        return flag_next

    def generate_consensus(self, th_merge=0.05, th_noisy=0.01, th_synchronize=0.3):   # TODO: reconsider the method
        '''Calculate a consensus sequence for each cluster.
        1. Initial consensus sequences are computed by Consed.
        2. Every two consensus sequences which have dis-similarity less than <th_merge> are merged.
           If you do not want to merge any clusters, then specify <th_merge=0>.
        3. Any cluster whose size is less than (<self.N> * <th_noisy>) is removed.
           Specify <th_noisy=0> to keep all clusters.
        4. Synchronize the phase of clusters which share dis-similarity less than <th_synchronize>.
        '''
        assert np.max(self.assignment) >= 0, 'No clustering result yet'

        # Initial consensus sequences
        self.cons_seqs = self._generate_consensus()
        logger.info(f'Unsynchronized consensus sequences:\n{self.cons_seqs}')

        # Merge too close clusters
        n_clusters_prev = self.cons_seqs.shape[0]
        while self._merge_clusters(th_merge):
            pass
        if self.cons_seqs.shape[0] != n_clusters_prev:
            logger.info(f'Merged similar clusters:\n{self.cons_seqs}')

        # Remove remaining noisy clusters
        del_row = [index for index, df in self.cons_seqs.iterrows()
                   if df['cluster_size'] < max([2, self.N * th_noisy])]   # too small cluster
        if len(del_row) > 0:
            self.cons_seqs = self.cons_seqs.drop(del_row).reset_index(drop=True)
            logger.debug(f'Removed too small clusters:\n{self.cons_seqs}')

        # Synchronize phase of the consensus units (= start position)
        n_cons = self.cons_seqs.shape[0]
        if n_cons == 1:
            return
        for i in range(n_cons - 1):   # TODO: simultaneously synchronize, or fix single seed
            for j in range(i + 1, n_cons):
                align = self.runner.align(self.cons_seqs['sequence'].iloc[i],
                                          self.cons_seqs['sequence'].iloc[j])
                if align.diff < th_synchronize:
                    logger.debug(f'Synchronize {i} and {j} (strand = {align.strand})')
                    self.cons_seqs.loc[j, 'sequence'] = align.mapped_seq(self.cons_seqs['sequence'].iloc[j])
        logger.info(f'Synchronized consensus sequences:\n{self.cons_seqs}')


if __name__ == '__main__':
    '''Only for internal usage by calc_dist_mat.'''
    p = argparse.ArgumentParser()
    p.add_argument('clustering_obj_pkl', type=str)
    p.add_argument('rows_pkl', type=str)
    p.add_argument('out_pkl', type=str)
    p.add_argument('n_core', type=int)
    args = p.parse_args()

    c = load_pickle(args.clustering_obj_pkl)
    save_pickle(c._calc_dist_mat(load_pickle(args.rows_pkl), args.n_core), args.out_pkl)
