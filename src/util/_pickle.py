from os.path import getsize
from typing import Any
import pickle
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
