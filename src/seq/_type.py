from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

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
class SeqInterval:
    """Class for an interval on a sequence."""

    b: int
    e: int

    @property
    def length(self) -> int:
        return self.e - self.b


@dataclass
class BedRecord(ExplicitRepr):
    """Class for an interval on a chromosomal sequence."""

    chr: str
    b: int
    e: int

    def __repr__(self) -> str:
        return self._order_repr(list(vars(self).keys()))

    @property
    def length(self) -> int:
        return self.e - self.b


@dataclass
class SatRecord(BedRecord):
    unit_seq: str
    n_copy: float

    @property
    def array_len(self):
        return self.e - self.b

    @property
    def unit_len(self):
        return len(self.unit_seq)
