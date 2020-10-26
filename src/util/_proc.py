from typing import Optional
import subprocess as sp
from multiprocessing import Process
from multiprocessing.pool import Pool
from logzero import logger


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
