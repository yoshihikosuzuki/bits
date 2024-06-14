from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List

import numpy as np


class ExplicitRepr(ABC):
    """An abstract class where the order of variables displayed in the repr()
    string must be explicityly specified.

    The `__repr__` method must be defined in every child class to specify the
    order of variables to be displayed in the repr() string.
    Pass a list of variable names to the helper method, `_order_repr()`, to
    generate the repr() string.
    """

    def _order_repr(self, var_names: List[str]) -> str:
        var_reprs = ", ".join(map(lambda x: f"{x}={repr(getattr(self, x))}", var_names))
        return f"{self.__class__.__name__}({var_reprs})"

    @abstractmethod
    def __repr__(self) -> str:
        return self._order_repr(["seq"])  # NOTE: This is an example


@dataclass
class SeqRecord(ExplicitRepr):
    """Abstract class for sequence object."""

    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
class FastaRecord(SeqRecord):
    """Sequence with name."""

    name: str

    def __repr__(self) -> str:
        return self._order_repr(["name", "seq"])


@dataclass
class FastqRecord(FastaRecord):
    """Sequence with name and base qualities."""

    qual: str

    def __repr__(self) -> str:
        return self._order_repr(["name", "seq", "qual"])

    @property
    def qual_phred(self) -> np.ndarray:
        return np.array(list(map(lambda c: ord(c) - 33, self.qual)), dtype=np.int8)


@dataclass
class DazzRecord(FastaRecord):
    """Sequence with name and DAZZ_DB ID."""

    id: int

    def __repr__(self) -> str:
        return self._order_repr(["id", "name", "seq"])


@dataclass
class SeqInterval(ExplicitRepr):
    """Class for an interval on a sequence."""

    b: int
    e: int

    @property
    def length(self) -> int:
        return self.e - self.b


@dataclass
class BedRecord(SeqInterval):
    """Class for an interval on a chromosomal sequence."""

    chr: str

    def __repr__(self) -> str:
        return self._order_repr(["chr", "b", "e"])

    @classmethod
    def from_string(cls, region: str):
        """Convert from e.g. `chr1:1000-2000` (1-index, closed) into 0-index, open"""
        chrom, b_e = region.split(":")
        b, e = b_e.split("-")
        return cls(chr=chrom, b=int(b) - 1, e=int(e))

    def to_string(self, comma: bool = False):
        """Covert 0-index, end open into 1-index, end closed (e.g. `chr1:101-200`)"""
        if not comma:
            return f"{self.chr}:{self.b + 1}-{self.e}"
        else:
            return f"{self.chr}:{self.b + 1:,}-{self.e:,}"

    @property
    def length(self) -> int:
        return self.e - self.b

    @property
    def string(self) -> str:
        return self.to_string()


@dataclass
class SatRecord(BedRecord):
    unit_seq: str
    n_copy: float

    def __repr__(self) -> str:
        return self._order_repr(["chr", "b", "e", "unit_seq", "n_copy"])

    @property
    def array_len(self):
        return self.length

    @property
    def unit_len(self):
        return len(self.unit_seq)
