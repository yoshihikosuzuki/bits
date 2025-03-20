"""
Defining two most basic classes, SeqRecord and SegRecord, and their inherited classes.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from logzero import logger

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
# Segment classes
##################################################################################################


# NOTE: implemented as a native class because of default values in SegRecord and
#       non-default values in BedRecord
class SegRecord:
    """A class for a segment/region object, which can be:
        - only a chromosome name (segment of the entire chromosome)
        - only begin/end positions (segment on an anonymous sequence)
        - both chromosome namd and begin/end positions (segment on the chromsome)

    Also, optionally the sequence information on which the segment is defined can be
    associated as `_seq`, whose type is either string or `SeqRecord`.

    `b` and `e` are pythonic (i.e. 0-indexed, start-closed, end-open).
    """

    def __init__(
        self,
        chr: Optional[str] = None,
        b: Optional[int] = None,
        e: Optional[int] = None,
        # _seq: Optional[Union[str, SeqRecord]] = None,
    ):
        self.chr = chr
        self.b = b
        self.e = e
        # self._seq = _seq

        if (self.b is None) != (self.e is None):
            logger.error("`b` and `e` must be both None or both non-None.")
        elif self.chr is None and self.b is None:
            logger.error("Either `chr` or `b` (and `e`) must be specified.")

    def __repr__(self) -> str:
        attr_names = ["chr", "b", "e"]
        text = ", ".join(
            [
                f"{attr_name}={repr(getattr(self, attr_name))}"
                for attr_name in attr_names
            ]
        )
        return f"{self.__class__.__name__}({text})"

    @classmethod
    def from_string(cls, region: str):
        """Convert from e.g. `chr1:100-200` (1-index, closed) into 0-index, open"""
        data = region.split(":")
        if len(data) == 1:
            # only chromosome, e.g. "chr1"
            chrom, b, e = data[0], None, None
        else:
            chrom, b_e = data
            data = b_e.split("-")
            if len(data) == 1:
                # single position is allowed only when converting from string
                # e.g. "chr1:100" is regarded as abbreviation of "chr1:100-100"
                b, e = int(data[0]), int(data[0])
            else:
                # e.g. "chr1:100-200"
                b, e = map(int, data)
            b -= 1
        return cls(chr=chrom, b=b, e=e)

    def to_string(self, comma: bool = False) -> str:
        """Covert 0-index, end open into 1-index, end closed (e.g. `chr1:100-200`)"""
        if self.chr is None:
            logger.error("Cannot convert to a string because `chr` is None.")
        elif self.b is None:
            return f"{self.chr}"
        elif not comma:
            return f"{self.chr}:{self.b + 1}-{self.e}"
        else:
            return f"{self.chr}:{self.b + 1:,}-{self.e:,}"

    @property
    def length(self) -> int:
        if self.b is None:
            if self._seq is not None:
                return self._seq.length
            else:
                logger.error("Entire chromosome, but no information about it.")
        else:
            return self.e - self.b

    # @property
    # def seq(self) -> str:
    #     if self._seq is None:
    #         logger.error("Sequence is not associated.")
    #     seq = self._seq if isinstance(self._seq, str) else self._seq.seq
    #     if self.b is None:
    #         # entire sequence
    #         return seq
    #     else:
    #         return seq[self.b : self.e]


@dataclass
class BedRecord(SegRecord):
    """Segment with a restriction that all `chr`, `b`, and `e` must be specified.
    Also, attributes can be specified.
    """

    chr: str
    b: int
    e: int

    def __init__(self, chr: str, b: int, e: int):
        super().__init__(chr, b, e)

    def __post_init__(self):
        if self.chr is None or self.b is None or self.e is None:
            logger.error("All of `chr`, `b`, and `e` must be specified.")

    # Modify __repr__ so all the attributes are displayed
    def __repr__(self) -> str:
        fixed_attr_names = ["chr", "b", "e"]
        attr_names = fixed_attr_names + list(
            set(vars(self).keys()) - set(fixed_attr_names)
        )
        text = ", ".join(
            [
                f"{attr_name}={repr(getattr(self, attr_name))}"
                for attr_name in attr_names
            ]
        )
        return f"{self.__class__.__name__}({text})"


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
class GffRecord(BedRecord):
    forward: bool
    type: str
    source: str

    # Modify __repr__ so all the attributes are displayed
    def __repr__(self) -> str:
        fixed_attr_names = ["chr", "b", "e", "forward", "type"]
        attr_names = fixed_attr_names + list(
            set(vars(self).keys()) - set(fixed_attr_names)
        )
        text = ", ".join(
            [
                f"{attr_name}={repr(getattr(self, attr_name))}"
                for attr_name in attr_names
            ]
        )
        return f"{self.__class__.__name__}({text})"
