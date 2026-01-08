from typing import List, Sequence, Tuple, Union

import numpy as np


def change_case(seq: str, case: str) -> str:
    assert case in (
        "original",
        "lower",
        "upper",
    ), "`case` must be 'original', 'lower', or 'upper'"
    return (
        seq if case == "original" else seq.lower() if case == "lower" else seq.upper()
    )


def findall(seq: str, query: str) -> List[int]:
    """List up all the positions of `query` in `seq`."""
    pos = []
    offset = 0
    while True:
        i = seq.find(query)
        if i == -1:
            break
        pos.append(offset + i)
        seq = seq[i + 1 :]
        offset += i + 1
    return pos


def split_at(
    seq: str,
    pos: Union[int, Sequence[int]],
    b: int = 0,
) -> List[str]:
    """Split a single string into a list of strings divided at specified position(s).

    @ seq
        Input string.
    @ pos
        Position(s) where `seq` is divided.
    @ b
        Use this option when `seq = _seq[b:b+len(seq)]` and `pos` is between `b` and `b+len(seq)`.
    """
    if isinstance(pos, int):
        pos = [pos]
    pos = [0] + [p - b for p in pos] + [len(seq)]
    return [seq[pos[i] : pos[i + 1]] for i in range(len(pos) - 1)]


def split_seq(seq: str, width: int) -> List[str]:
    """Split `seq` at every `width` characters."""
    return [seq[i : i + width] for i in range(0, len(seq), width)]


def reverse_seq(seq: str) -> str:
    return seq[::-1]


F_TO_RC = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def revcomp_seq(seq: str) -> str:
    return "".join([F_TO_RC[c] for c in seq[::-1]])


def run_length_encoding(seq: str) -> List[Tuple[str, int]]:
    if len(seq) == 0:
        return []

    encoded = []
    prev_char = seq[0]
    count = 1
    for char in seq[1:]:
        if char == prev_char:
            count += 1
        else:
            encoded.append((prev_char, count))
            prev_char = char
            count = 1
    encoded.append((prev_char, count))

    return encoded


def compress_homopolymer(seq: str) -> str:
    comp_seq, prev_base = "", ""
    for base in seq:
        if base != prev_base:
            comp_seq += base
        prev_base = base
    return comp_seq


def ascii_to_phred(c: str) -> int:
    """Convert quality character in fastq into Phred score."""
    assert len(c) == 1, "`c` must be a single character"
    return ord(c) - 33


def phred_to_log10_p_error(phred: int) -> float:
    """Convert Phred score into log10(Pr{base is erroneous})."""
    assert 0 <= phred <= 93, "`phred` must be in a range of [0..93]"
    return -phred / 10


PHRED_TO_LOG_CORRECT = [
    -np.inf if phred == 0 else np.log10(1 - np.power(10, -phred / 10))
    for phred in range(94)
]


def phred_to_log10_p_correct(phred: int) -> float:
    """Convert Phred score into log10(Pr{base is correct})."""
    assert 0 <= phred <= 93, "`phred` must be in a range of [0..93]"
    return PHRED_TO_LOG_CORRECT[phred]
