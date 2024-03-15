import subprocess as sp
import sys
from multiprocessing import Process
from multiprocessing.context import DefaultContext
from multiprocessing.pool import Pool

from logzero import logger


def run_command(command: str, err_msg: str = "", exit_on_error: bool = True) -> str:
    """General-purpose shell command runner.

    positional arguments:
      @ command: A string evaluated as a line in bash.
      @ err_msg: Any text shown on error.
      @ exit_on_error: If True, exit the process on error.
    """
    try:
        out = sp.check_output(command, shell=True)
    except sp.CalledProcessError as proc:
        if err_msg != "":
            err_msg = f"({err_msg})"
        logger.exception(proc)
        logger.error(f"Command failed: {command} {err_msg}")
        if exit_on_error:
            sys.exit(1)
        return proc.output.decode("utf-8")
    else:
        return out.decode("utf-8")


class NoDaemonPool(Pool):
    """Custom Pool that can have child processes."""

    def get_context(self):
        return NoDaemonContext()


class NoDaemonContext(DefaultContext):
    def Process(self, *args, **kwargs):
        return NoDaemonProcess(*args, **kwargs)


class NoDaemonProcess(Process):
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)
