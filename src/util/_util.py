from os.path import getsize
from typing import Any, Optional
import pickle
import subprocess as sp
from multiprocessing import Process
from multiprocessing.pool import Pool
from logzero import logger


def save_pickle(obj: Any,
                out_pkl_fname: str,
                protocol: int = 4,
                verbose: bool = True) -> None:
    """Utility for pickling and saving an object.

    positional arguments:
      @ obj           : Serializable object.
      @ out_pkl_fname : Output file name.

    optional arguments:
      @ protocol : For pickling.
      @ verbose  : If True, show logs.
    """
    ret = pickle.dumps(obj, protocol=protocol)
    with open(out_pkl_fname, 'wb') as f:
        f.write(ret)
    if verbose:
        length = f"(length = {len(obj)}) " if hasattr(obj, '__len__') else ""
        logger.info(f"Saved {type(obj)} object {length}"
                    f"to {out_pkl_fname} ({len(ret)} bytes)")


def load_pickle(pkl_fname: str,
                verbose: bool = True) -> Any:
    """Utility for loading and unpickling an object.

    positional arguments:
      @ pkl_fname : Input file name.

    optional arguments:
      @ verbose  : If True, show logs.
    """
    with open(pkl_fname, 'rb') as f:
        ret = pickle.load(f)
    if verbose:
        length = f"(length = {len(ret)}) " if hasattr(ret, '__len__') else ""
        logger.info(f"Loaded {type(ret)} object {length}"
                    f"from {pkl_fname} ({getsize(pkl_fname)} bytes)")
    return ret


def run_command(command: str,
                err_msg: Optional[str] = None) -> str:
    """General-purpose shell command runner.

    positional arguments:
      @ command : A string evaluated as a line in bash.
    """
    try:
        out = sp.check_output(command, shell=True)
    except sp.CalledProcessError as proc:
        _err_msg = f"({err_msg})" if err_msg is not None else ""
        logger.error(f"Error raised with command: {command} {_err_msg}")
        logger.exception(proc)
        return proc.output.decode("utf-8")
    else:
        return out.decode("utf-8")


class NoDaemonProcess(Process):
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class NoDaemonPool(Pool):
    """Custom Pool that can have child processes."""
    Process = NoDaemonProcess
