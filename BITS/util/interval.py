# Utilities for pyinterval package
from interval import interval


def interval_len(intvls):
    """
    Calculate the sum of the length of the intervals in <intvls>.
    NOTE: It is assumed that <intvls> is INTEGER interval object.
    """

    return sum([i[0][1] - i[0][0] + 1 for i in intvls.components])


def subtract_interval(a, b, th_length=0):
    """
    Calculate A - B as interval objects and then remove remaining intervals shorter than <th_length>.
    NOTE: It is assumed that <a> and <b> are INTEGER interval objects.
    """
    
    ret = interval()
    ai, bi = 0, 0
    an, bn = len(a), len(b)

    while ai < an and bi < bn:
        al, ar = a[ai]
        while bi < bn and b[bi][1] < al:   # skip non-involved intervals
            bi += 1
        while bi < bn and b[bi][0] < ar:
            bl, br = b[bi]
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
        al, ar = a[ai]
        ret |= interval([al, ar])
        ai += 1

    # Remove intervals shorter than <th_length>
    return interval(*[(i[0][0], i[0][1]) for i in ret.components if interval_len(i) >= th_length])
