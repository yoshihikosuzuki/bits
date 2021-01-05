from typing import List
import numpy as np


def findall(seq: str, query: str) -> List[int]:
    pos = []
    offset = 0
    while True:
        i = seq.find(query)
        if i == -1:
            break
        pos.append(offset + i)
        seq = seq[i + 1:]
        offset += i + 1
    return pos


def split_seq(seq: str, width: int) -> List[str]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def reverse_seq(seq: str) -> str:
    return seq[::-1]


F_TO_RC = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def revcomp_seq(seq: str) -> str:
    return ''.join([F_TO_RC[c] for c in seq[::-1]])


def compress_homopolymer(seq: str) -> str:
    comp_seq, prev_base = "", ''
    for base in seq:
        if base != prev_base:
            comp_seq += base
        prev_base = base
    return comp_seq


def ascii_to_phred(c: str) -> int:
    assert len(c) == 1, "`c` must be a single character"
    return ord(c) - 33


def phred_to_log10_p_error(phred: int) -> float:
    assert 0 <= phred <= 93, "`phred` must be in a range of [0..93]"
    return -phred / 10


PHRED_TO_LOG_CORRECT = [-np.inf if phred == 0
                        else np.log10(1 - np.power(10, -phred / 10))
                        for phred in range(94)]


def phred_to_log10_p_correct(phred: int) -> float:
    assert 0 <= phred <= 93, "`phred` must be in a range of [0..93]"
    return PHRED_TO_LOG_CORRECT[phred]
