from dataclasses import dataclass
from os.path import isfile
from typing import List, Optional, Sequence, Union

import bits.util as bu

from ._type import SeqRecord


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


def calc_lens(
    data: Union[str, Sequence[Union[SeqRecord, str, int]]], force: bool = False
) -> List[int]:
    """Utility function for obtaining sequence lengths.

    Parameters
    ----------
    data
        Input data. Must be one of the followings:
            1) fasta/q file name
            2) list of sequence objects
            3) list of sequences
            4) list of sequence lengths

    Returns
    -------
        List of sequence lengths.
    """
    lens = None
    if isinstance(data, str):  # Fastx file
        in_fastx = data
        out_nl = f"{in_fastx}.nl"
        if not isfile(out_nl) or force:
            bu.run_command(f"seqkit fx2tab -nl {in_fastx} >{out_nl}")
        lens = list(
            map(
                int,
                bu.run_command(f"cut -f2 {out_nl}").strip().split("\n"),
            )
        )
    elif isinstance(data, Sequence):
        if isinstance(data[0], SeqRecord):  # list of sequence objects
            seqs = data
            lens = [seq.length for seq in seqs]
        elif isinstance(data[0], str):  # list of sequences
            seqs = data
            lens = [len(seq) for seq in seqs]
        elif isinstance(data[0], int):  # list of sequence lengths
            lens = data
    assert lens is not None, "Failed to guess the type of input data"
    return lens


def calc_nx(
    data: Union[str, Sequence[Union[SeqRecord, str, int]]],
    G: Optional[int] = None,
    force: bool = False,
) -> List[int]:
    """Calaulate Nx (or NGx if G is given) values from a fasta file.

    Parameters
    ----------
    data
        Input data. Must be one of the followings:
            1) fasta/q file name
            2) list of sequence objects
            3) list of sequences
            4) list of sequence lengths
    G, optional
        Genome size, by default None

    Returns
    -------
        Nx/NGx values.
    """
    lens = calc_lens(data, force)
    if G is None:
        G = sum(lens)
    thres = [G * x / 100 for x in range(100 + 1)]
    thres_idx = 0
    nx = [0] * len(thres)
    s = 0
    for l in sorted(lens, reverse=True):
        s += l
        while thres_idx <= 100 and s >= thres[thres_idx]:
            nx[thres_idx] = l
            thres_idx += 1
    return nx
