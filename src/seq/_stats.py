from dataclasses import dataclass
from os.path import isfile
from typing import List, Optional

import bits.util as bu


@dataclass
class SeqStats:
    fname: str
    nseqs: int
    nbases: int
    min_len: int
    mean_len: int
    max_len: int
    n50_len: int


def load_seq_stats(in_fname: str) -> List[SeqStats]:
    assert isfile(in_fname), f"Input file ({in_fname}) not found"
    ret = []
    with open(in_fname) as f:
        f.readline()
        for line in f:
            data = line.strip().split()
            fname = data[0]
            nseqs = int(data[3].replace(",", ""))
            nbases = int(data[4].replace(",", ""))
            min_len = int(data[5].replace(",", ""))
            mean_len = round(float(data[6].replace(",", "")))
            max_len = int(data[7].replace(",", ""))
            n50_len = int(data[12].replace(",", "")) if len(data) > 12 else None
            ret.append(
                SeqStats(fname, nseqs, nbases, min_len, mean_len, max_len, n50_len)
            )
    return ret


def calc_seq_stats(in_fastx: str, force: bool = False) -> List[SeqStats]:
    out_stats = f"{in_fastx}.stats"
    if not isfile(out_stats) or force:
        bu.run_command(f"seqkit stats -a {in_fastx} >{out_stats}")
    return load_seq_stats(out_stats)


def calc_nx(in_fastx: str, G: Optional[int] = None) -> List[int]:
    """Calaulate Nx (or NGx if G is given) values from a fasta file.

    Args:
        in_fastx (str): Input fasta file.
        G (Optional[int], optional): Genome size. Defaults to None.

    Returns:
        List[int]: Nx/NGx values.
    """
    lens = list(
        map(
            int,
            bu.run_command(f"seqkit fx2tab -nl {in_fastx} | cut -f2")
            .strip()
            .split("\n"),
        )
    )
    if G is None:
        G = sum(lens)
    thres = [G * x / 100 for x in range(100 + 1)]
    thres_idx = 0
    nx = [0] * len(thres)
    s = 0
    for l in sorted(lens, reverse=True):
        s += l
        while s > thres[thres_idx]:
            nx[thres_idx] = l
            thres_idx += 1
    return nx
