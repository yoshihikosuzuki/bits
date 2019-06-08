from functools import wraps
import logging
import logzero
from logzero import logger


def suppress_debug_log():
    logzero.loglevel(logging.INFO)


def print_log(proc_name, show_args=True):
    """
    Simple decolator for watching start and end of a function.
    """
    def _print_log(func):
        @wraps(func)   # necessary to show information of the original function when doing `func?`
        def wrapper(*args, **kwds):
            args_info = f"(args={args}, kwds={kwds})" if show_args else ""
            logger.info(f"Starting {proc_name} {args_info}")
            ret = func(*args, **kwds)
            logger.info(f"Finished {proc_name}")
            return ret
        return wrapper
    return _print_log
