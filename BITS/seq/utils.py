RC_MAP = dict(zip("ACGTacgtNn-", "TGCAtgcaNn-"))


def revcomp(seq):
    return ''.join([RC_MAP[c] for c in seq[::-1]])


def reverse_seq(seq):
    return ''.join(list(reversed(seq)))


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

