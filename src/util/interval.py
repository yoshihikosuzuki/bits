"""Utilities for pyinterval package [https://github.com/taschini/pyinterval]."""
from interval import interval


def intvl_len(intvls):
    """Calculate the sum of the interval lengths of <intvls>.
    <intvls> must consist of INTEGER intervals with inclusive boundaries.
    """
    return sum([i[0][1] - i[0][0] + 1 for i in intvls.components])


def subtract_intvl(a_intvl, b_intvl, th_length=0):
    """Calculate <a_intvl> - <b_intvl> as interval.
    <intvls> must consist of INTEGER intervals with inclusive boundaries.
    """
    ret = interval()
    ai, bi, an, bn = 0, 0, len(a_intvl), len(b_intvl)
    while ai < an and bi < bn:
        al, ar = a_intvl[ai]
        while bi < bn and b_intvl[bi][1] < al:   # skip non-involved intervals
            bi += 1
        while bi < bn and b_intvl[bi][0] < ar:
            bl, br = b_intvl[bi]
            if ar < bl:   # no overlap
                break
            if al < bl:
                ret |= interval([al, bl - 1])
            if br < ar:   # interval remains
                al = br + 1   # shorten the interval
            else:
                al = ar + 1
                break
            bi += 1
        if al <= ar:
            ret |= interval([al, ar])
        ai += 1
    while ai < an:   # remaining intervals
        al, ar = a_intvl[ai]
        ret |= interval([al, ar])
        ai += 1
    return ret


def filter_short_intvls(intvls, min_len):
    """Filter out intervals in <intvls> shorter than <min_len>.
    <intvls> must consist of INTEGER intervals with inclusive boundaries.
    """
    return interval(*[i for i in intvls if intvl_len(i) >= min_len])
