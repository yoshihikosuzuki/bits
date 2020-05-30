import subprocess as sp
from multiprocessing import Process
from multiprocessing.pool import Pool
from logzero import logger


def run_command(command: str) -> str:
    """General-purpose shell command executer."""
    try:
        out = sp.check_output(command, shell=True)
    except sp.CalledProcessError as proc:
        logger.warn(f"Error raised during: {command}")
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
