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
             out_log="sge_stdout",
             err_log="sge_stderr",
             n_core=1,
             time_limit=None,   # NOTE: not supported yet
             mem_limit=None,   # NOTE: not supported yet
             depend=[],
             wait=False):
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
    if depend is not None:
        header += f"\n#$ -hold_jid {depend}"

    return f"{header}\n\n{script}\n"


def submit_sge(script,
               out_fname,
               job_name="run_script",
               out_log="sge.log",
               n_core=1,
               depend=[],
               wait=False,
               submit_command="qsub"):
    """
    Submit a script as a SGE job with <submit_command> after adding headers with specified options.
    """

    with open(out_fname, 'w') as f:
        f.write(sge_nize(script,
                         job_name,
                         out_log,
                         n_core,
                         depend,
                         wait))

    return run_command(f"{submit_command} {out_fname}").split()[2]   # job ID


def slurm_nize(script,
               job_name="run_script",
               out_log="sbatch_stdout",
               err_log="sbatch_stderr",
               n_core=1,
               time_limit="24:00:00",
               mem_per_cpu=1024,
               partition="batch",
               depend=None,
               wait=False):
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
    if depend is not None:
        header += f"\n#SBATCH -d {depend}"

    return f"{header}\n\n{script}\n"


def submit_slurm(script,
                 job_name="run_script",
                 out_log="sbatch_stdout",
                 err_log="sbatch_stderr",
                 n_core=1,
                 time_limit="24:00:00",
                 mem_per_cpu=1024,
                 partition="batch",
                 depend=None,
                 wait=False):
    """
    Add headers for sbatch of SLURM.
    """


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


def subtract_interval(a, b, th_length=0):
    """
    For pyinterval (https://pyinterval.readthedocs.io/en/latest/).
    Calculate A - B, where A and B are interval objects.
    Afterwards, remove remaining intervals shorter than <th_length>.
    """
    
    from interval import interval

    ret = interval()
    ai, bi = 0, 0
    an, bn = len(a), len(b)

    while ai < an and bi < bn:
        al, ar = a[ai]
        while bi < bn and b[bi][1] < al:   # skip non-involved intervals
            bi += 1
        while bi < bn and b[bi][0] < ar:
            bl, br = b[bi]
            if ar < bl:   # no overlap
                break
            if al < bl:
                ret |= interval([al, bl - 1])
            if br < ar:   # interval remains
                al = br + 1   # shorten the interval
            else:
                al = ar + 1
                break
            bi += 1
        if al <= ar:
            ret |= interval([al, ar])
        ai += 1

    while ai < an:   # remaining intervals
        al, ar = a[ai]   # a[ai] is Component, not interval
        ret |= interval([al, ar])
        ai += 1

    if th_length == 0:
        return ret
    else:
        # Remove intervals shorter than <th_length>
        return interval(*[(i[0][0], i[0][1]) for i in ret.components if interval_len(i) >= th_length])
