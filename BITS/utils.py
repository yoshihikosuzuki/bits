from multiprocessing import Process
from multiprocessing.pool import Pool
from logzero import logger


def load_pickle(pkl_fname):
    import pickle

    with open(pkl_fname, 'rb') as f:
        return pickle.load(f)


def save_pickle(obj, pkl_fname):
    import pickle

    with open(pkl_fname, 'wb') as f:
        pickle.dump(obj, f)


def run_command(command, show_error_msg=False):
    """
    General-purpose command executer.
    """

    import subprocess

    try:
        out = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as proc:
        logger.warn(f"An error occured while command execution")
        if show_error_msg:
            logger.exception(proc)
        return proc.output.decode('utf-8')
    else:
        return out.decode('utf-8')


def print_log(process_name, show_args=True):
    """
    Simple decolator for watching start and end of a function.
    """

    def _print_log(func):
        def wrapper(*args, **kwds):
            args_info = f"(args={args}, kwds={kwds})" if show_args else ""
            logger.info(f"Starting {process_name} {args_info}")
            ret = func(*args, **kwds)
            logger.info(f"Finished {process_name}")
            return ret
        return wrapper
    return _print_log


def sge_nize(script,
             job_name="run_script",
             out_log="sge.log",
             n_core=1,
             wait=True):
    """
    Add headers for qsub of SGE.
    """

    header = '\n'.join([f"#!/bin/bash",
                        f"#$ -N {job_name}",
                        f"#$ -o {out_log}",
                        f"#$ -j y",
                        f"#$ -S /bin/bash",
                        f"#$ -cwd",
                        f"#$ -V",
                        f"#$ -pe smp {n_core}",
                        f"#$ -sync {'y' if wait is True else 'n'}"])

    return f"{header}\n\n{script}\n"


def slurm_nize(script,
               job_name="run_script",
               out_log="sbatch_stdout",
               err_log="sbatch_stderr",
               n_core=1,
               time_limit="24:00:00",
               mem_per_cpu=1024,
               partition="batch",
               wait=True):
    """
    Add headers for sbatch of SLURM.
    """

    header = '\n'.join([f"#!/bin/bash",
                        f"#SBATCH -J {job_name}",
                        f"#SBATCH -o {out_log}",
                        f"#SBATCH -e {err_log}",
                        f"#SBATCH -n 1",
                        f"#SBATCH -N 1",
                        f"#SBATCH -c {n_core}",
                        f"#SBATCH -t {time_limit}",
                        f"#SBATCH --mem-per-cpu={mem_per_cpu}",
                        f"#SBATCH --partition={partition}",
                        f"{'#SBATCH --wait' if wait is True else ''}"])

    return f"{header}\n\n{script}\n"


class NoDaemonProcess(Process):
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class NoDaemonPool(Pool):
    """
    Inherited class of Pool so that it runs as a non-daemon process and thus can have child processes.
    """
    Process = NoDaemonProcess


def make_line(x0, y0, x1, y1, col='black', width=1, layer='below'):
    """
    For Plotly.
    Create a line-shape object for Plotly.
    """

    return {'type': 'line', 'xref': 'x', 'yref': 'y',
            'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
            'line': {'color': col, 'width': width},
            'layer': layer}


def interval_len(intvls):
    """
    For pyinterval (https://pyinterval.readthedocs.io/en/latest/).
    Return the sum of the interval lengths in <intvls>.
    """
    
    ret = 0
    for intvl in intvls.components:
        ret += intvl[0][1] - intvl[0][0] + 1
    return ret


def subtract_interval(a_interval, b_interval, length_threshold=0):
    """
    For pyinterval (https://pyinterval.readthedocs.io/en/latest/).
    Calculate A - B, where A and B are interval objects.
    Afterwards, remove remaining intervals shorter than <length_threshold>.
    """
    
    from interval import interval

    ret_interval = interval()
    intersection = a_interval & b_interval
    a_index = 0
    i_index = 0

    flag_load_new_interval = True
    while a_index < len(a_interval) and i_index < len(intersection):
        if flag_load_new_interval:
            a_intvl = a_interval[a_index]
        else:
            flag_load_new_interval = False

        if a_intvl[1] < intersection[i_index][0]:
            ret_interval |= interval(a_intvl)
            a_index += 1
        elif a_intvl[0] > intersection[i_index][1]:
            i_index += 1
        else:
            a_start, a_end = a_intvl
            i_start, i_end = intersection[i_index]
            start_contained = True if min(a_start, i_start) == i_start else False
            end_contained = True if max(a_end, i_end) == i_end else False
            if start_contained:
                if end_contained:
                    a_index += 1
                    i_index += 1
                else:
                    # the tail interval that did not intersect
                    # with this B-interval (but maybe with the next B-interval)
                    a_intvl = interval[i_end + 1, a_end]
                    flag_load_new_interval = False
                    i_index += 1
            else:
                if end_contained:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_index += 1
                else:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_intvl = interval[i_end + 1, a_end]
                    flag_load_new_interval = False
                    i_index += 1
    if not flag_load_new_interval:   # add a tail interval if it exists
        ret_interval |= interval(a_intvl)
        a_index += 1
    while a_index < len(a_interval):   # add remaining intervals in A
        ret_interval |= interval(a_interval[a_index])
        a_index += 1

    ret_interval_long = interval()   # length threshold
    for intvl in ret_interval.components:
        if interval_len(intvl) >= length_threshold:
            ret_interval_long |= intvl

    return ret_interval_long
