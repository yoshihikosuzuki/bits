from __future__ import annotations
from typing import Iterator, Tuple

CIGAR_CHARS = {'=', 'X', 'I', 'D'}


class Cigar(str):
    """A class for a CIGAR string that can be treated as a str.

    usage:
      > cigar = Cigar("15=1X2D3=")
      > cigar
      '15=1X2D3='
      > list(cigar)
      [(15, '='), (1, 'X'), (2, 'D'), (3, '=')]
      > for l, op in cigar:
      >     print(l, op)
      15 =
      1 X
      2 D
      3 =

    positional arguments:
      @ cigar <str> : CIGAR string
    """
    def __iter__(self) -> Iterator[Tuple[int, str]]:
        """Iterator over tuples of length and operation."""
        length = ""
        for c in super().__iter__():
            # if c in CIGAR_CHARS:
            if not c.isnumeric():
                yield (int(length), c) if length != "" else (1, c)
                length = ""
            else:
                length += c

    def __reversed__(self) -> Iterator[Tuple[int, str]]:
        """Used for reverse or reverse complement."""
        yield from reversed(list(self))

    @property
    def aln_length(self) -> int:
        """Alignment length including gaps."""
        return sum([l for l, _ in self])

    def reverse(self) -> Cigar:
        """Reverse CIGAR without swapping I/D. Used for reverse complement."""
        return Cigar(''.join(map(lambda x: ''.join(map(str, x)), reversed(self))))

    def revcomp(self) -> Cigar:
        return self.reverse()

    def swap(self) -> Cigar:
        """Swap I and D. This swaps the role of query and target."""
        return Cigar(self.replace('I', '?').replace('D', 'I').replace('?', 'D'))

    def flatten(self) -> FlattenCigar:
        """Convert to a flatten CIGAR, a sequence of each operation."""
        return FlattenCigar(''.join([op * l for l, op in self]))


class FlattenCigar(str):
    """A class for a flatten CIGAR string that can be treated as a str.
    This is easier to handle than Cigar.

    usage:
      > cigar = Cigar("15=1X2D3=")
      > fcigar = cigar.flatten()
      > fcigar
      '===============XDD==='
      > fcigar.unflatten()
      '15=1X2D3='

    positional arguments:
      @ fcigar <str> : A flatten CIGAR string
    """
    @property
    def aln_length(self) -> int:
        """Implemented just for consistency with Cigar."""
        return len(self)

    def reverse(self) -> FlattenCigar:
        """Reverse CIGAR without swapping I/D. Used for reverse complement."""
        return FlattenCigar(reversed(self))

    def revcomp(self) -> FlattenCigar:
        return self.reverse()

    def swap(self) -> FlattenCigar:
        """Swap I and D. This swaps the role of query and target."""
        return FlattenCigar(self.replace('I', '?').replace('D', 'I').replace('?', 'D'))

    def unflatten(self) -> Cigar:
        """Convert to a normal CIGAR string."""
        cigar = ""
        length = 0
        prev_c = self[0]
        for c in self:
            if c == prev_c:
                length += 1
            else:
                cigar += f"{length}{prev_c}"
                length = 1
                prev_c = c
        cigar += f"{length}{prev_c}"
        return Cigar(cigar)
