from logzero import logger


def debug_mode(yes):
    import logging
    import logzero
    if yes:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)


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
