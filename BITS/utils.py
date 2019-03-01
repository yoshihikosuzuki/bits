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
        pickle.dump(obj, f, protocol=4)


def run_command(command, show_error_msg=True):
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
