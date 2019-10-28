import numpy as np

RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def reverse_seq(seq):
    return seq[::-1]


def revcomp_seq(seq):
    return ''.join([RC_MAP[c] for c in seq[::-1]])


def compress_homopolymer(seq):
    ret = ""
    prev = ""
    for s in seq:
        if s != prev:
            ret += s
        prev = s
    return ret


def split_seq(seq, width):
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def asciis_to_phreds(qual):
    return np.array([ascii_to_phred(a) for a in qual], dtype=np.int8)


def ascii_to_phred(a):
    return ord(a) - 33


def phred_to_log10_p_error(phred):
    return -phred / 10


def phred_to_log10_p_correct(phred):
    return np.log10(1 - np.power(10, -phred / 10))
