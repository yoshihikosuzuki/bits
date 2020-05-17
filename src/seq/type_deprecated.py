from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Overlap:
    """An overlap between two reads.
    `a_read[a_start:a_end]` matchs `strand(b_read)[b_start:b_end]`.
    `[a|b]_read` can be any identifier of a read: ID/name/object/sequence/etc.
    """
    a_read: Any
    b_read: Any
    strand: int
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    length: int
    diff: float
    cigar: Optional[Cigar]

    def __post_init__(self):
        assert self.strand in (0, 1), "`strand` must be 0/1"
