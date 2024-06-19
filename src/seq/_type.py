from dataclasses import dataclass
from typing import Optional

import numpy as np

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
class FastaRecord(SeqRecord):
    """Sequence with name."""

    name: str
    seq: str = SeqRecord.__dataclass_fields__["seq"].default

    def __init__(self, name: str, seq: str):
        super().__init__(seq)
        self.name = name


@dataclass
class FastqRecord(FastaRecord):
    """Sequence with name and base qualities."""

    qual: str

    @property
    def qual_phred(self) -> np.ndarray:
        return np.array(list(map(lambda c: ord(c) - 33, self.qual)), dtype=np.int8)


@dataclass
class DazzRecord(FastaRecord):
    """Sequence with name and DAZZ_DB ID."""

    id: int
    name: str = FastaRecord.__dataclass_fields__["name"].default
    seq: str = FastaRecord.__dataclass_fields__["seq"].default

    def __init__(self, id: int, name: str, seq: str):
        super().__init__(name, seq)
        self.id = id


##################################################################################################
# Interval classes
##################################################################################################


@dataclass
class SeqInterval:
    """Abstract class for an interval object."""

    b: int
    e: int

    @property
    def length(self) -> int:
        return self.e - self.b


@dataclass
class BedRecord(SeqInterval):
    """Class for an interval on a chromosomal sequence."""

    chr: str
    b: int = SeqInterval.__dataclass_fields__["b"].default
    e: int = SeqInterval.__dataclass_fields__["e"].default

    def __init__(self, chr: str, b: int, e: int):
        super().__init__(b, e)
        self.chr = chr

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


# TODO: integrate into BedRecord?
@dataclass
class Region:
    """Class for a genomic region object.
    In contrast to `BedRecord`, `b` and `e` can be None."""

    chr: str
    b: Optional[int] = None
    e: Optional[int] = None

    @classmethod
    def from_string(cls, region: str):
        """Convert from e.g. `chr1:1000-2000` (1-index, closed) into 0-index, open"""
        data = region.split(":")
        if len(data) == 1:  # only chromosome. e.g. "chr1"
            chrom, b, e = data[0], None, None
        else:
            chrom, b_e = data
            data = b_e.split("-")
            if len(data) == 1:  # single position
                b, e = int(data[0]), int(data[0])
            else:
                b, e = map(int, data)
            b -= 1
        return cls(chr=chrom, b=b, e=e)


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
