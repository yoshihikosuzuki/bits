from typing import Callable, List

from ._proc import NoDaemonPool


def _run_parallel(func: Callable, args_list: List, multi_arg: bool = False) -> List:
    return [func(*args) if multi_arg else func(args) for args in args_list]


def run_parallel(
    func: Callable, args_list: List, n_proc: int, multi_arg: bool = False
) -> List:
    n_unit = -(-len(args_list) // n_proc)
    args_split = [
        (func, args_list[n_unit * i : n_unit * (i + 1)], multi_arg)
        for i in range(-(-len(args_list) // n_unit))
    ]
    with NoDaemonPool(n_proc) as pool:
        results = []
        for ret in pool.starmap(_run_parallel, args_split):
            results += ret
        return results
