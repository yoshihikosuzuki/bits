"""
Defining two most basic classes, SeqRecord and SegRecord, and their inherited classes.
"""

from dataclasses import dataclass
from typing import Optional, Union

import numpy as np
from logzero import logger

from ..util import ordered_repr, verbose_repr

##################################################################################################
# Sequence classes
##################################################################################################


@dataclass
class SeqRecord:
    """Abstract class for a sequence object."""

    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
@ordered_repr("name", "seq")
class FastaRecord(SeqRecord):
    """Sequence with name."""

    name: str


@dataclass
@ordered_repr("name", "seq", "qual")
class FastqRecord(FastaRecord):
    """Sequence with name and base qualities."""

    qual: str

    @property
    def qual_phred(self) -> np.ndarray:
        return np.array(list(map(lambda c: ord(c) - 33, self.qual)), dtype=np.int8)


@dataclass
@ordered_repr("id", "name", "seq")
class DazzRecord(FastaRecord):
    """Sequence with name and DAZZ_DB ID."""

    id: int


##################################################################################################
# Segment classes
##################################################################################################


# NOTE: SegRecord cannot be a dataclass because BedRecord inherits from it and
# overriding optional fields with non-optional fields, which is not allowed in
# dataclasses.
class SegRecord:
    """Abstract class for a segment/region/interval object.

    `b` and `e` are pythonic (0-indexed, start-closed, end-open).
    """

    def __init__(
        self,
        chrom: Optional[str] = None,
        b: Optional[int] = None,
        e: Optional[int] = None,
    ):
        self.chrom = chrom
        self.b = b
        self.e = e

    @classmethod
    def from_string(cls, region: str):
        """Convert from e.g. `chr1:100-200` (1-index, closed) into 0-index, open"""
        data = region.split(":")
        if len(data) == 1:
            chrom, b, e = data[0], None, None
        else:
            chrom, b_e = data
            data = b_e.split("-")
            if len(data) == 1:
                b, e = int(data[0]), int(data[0])
            else:
                b, e = map(int, data)
            b -= 1
        return cls(chrom=chrom, b=b, e=e)

    def to_string(self, comma: bool = False) -> str:
        """Convert 0-index,end-open into 1-index,end-closed (e.g. `chr1:100-200`)"""
        assert (
            self.chrom is not None
        ), "Cannot convert to a string because `chrom` is None."
        if self.b is None:
            return f"{self.chrom}"
        elif not comma:
            return f"{self.chrom}:{self.b + 1}-{self.e}"
        else:
            return f"{self.chrom}:{self.b + 1:,}-{self.e:,}"

    @property
    def length(self) -> int:
        assert self.b is not None and self.e is not None, logger.error(
            f"b and/or e is undefined"
        )
        return self.e - self.b


@dataclass
@verbose_repr("chrom", "b", "e")
class BedRecord(SegRecord):
    """Segment with a restriction that all `chrom`, `b`, and `e` must be specified."""

    chrom: str
    b: int
    e: int

    def __post_init__(self):
        assert (
            self.chrom is not None and self.b is not None and self.e is not None
        ), "All of `chrom`, `b`, and `e` must be specified."


@dataclass
class SatRecord(BedRecord):
    unit_seq: str
    n_copy: float

    @property
    def array_len(self):
        return self.length

    @property
    def unit_len(self):
        return len(self.unit_seq)


@dataclass
@verbose_repr("chrom", "b", "e", "forward", "type")
class GffRecord(BedRecord):
    forward: bool
    type: str
    source: str
