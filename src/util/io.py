from typing import Any
import pickle
from logzero import logger


def load_pickle(pkl_fname: str) -> Any:
    with open(pkl_fname, 'rb') as f:
        ret = pickle.load(f)
    length = f"(length = {len(ret)})" if hasattr(ret, '__len__') else ""
    logger.info(f"{pkl_fname}: Loaded {type(ret)} object {length}")
    return ret


def save_pickle(obj: Any, out_pkl_fname: str):
    with open(out_pkl_fname, 'wb') as f:
        pickle.dump(obj, f, protocol=4)
