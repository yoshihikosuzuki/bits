from dataclasses import dataclass
from typing import Any, Optional
import edlib
from .util import revcomp_seq, split_seq
from .cigar import Cigar


@dataclass(frozen=True)
class Alignment:
    """An overlap between two reads.

    optional arguments:
      @ [a|b]_seq          : Not limited to str but any identifier of a sequence.
                             e.g. ID/name/object/sequence/etc.
      @ strand             : 0 or 1
      @ apply_strand_first : `a_seq[a_start:a_end]` matches to:
                             `strand(b_seq)[b_start:b_end]` if True (default),
                             `strand(b_seq[b_start:b_end])` otherwise.
      @ [a|b]_[start|end]  : Start/end position of the alignment.
      @ [a|b]_len          : Length of the entire sequence.
      @ length             : Length of the alignment inlcuding gaps.
      @ diff               : Sequence dissimilarity of the alignment.
      @ cigar              : Cigar string of the alignment.
    """
    a_seq: Optional[Any] = None
    b_seq: Optional[Any] = None
    strand: int = 0
    apply_strand_first: bool = True
    a_start: Optional[int] = None
    a_end: Optional[int] = None
    a_len: Optional[int] = None
    b_start: Optional[int] = None
    b_end: Optional[int] = None
    b_len: Optional[int] = None
    length: Optional[int] = None
    diff: Optional[float] = None
    cigar: Optional[Cigar] = None

    def __post_init__(self):
        assert self.strand in (0, 1), "`strand` must be 0/1"

    def __repr__(self) -> str:
        a_repr = f"a[{self.a_start}:{self.a_end}]"
        b_range = f"[{self.b_start}:{self.b_end}]"
        b_repr = (f"b{b_range}" if self.strand == 0
                  else f"(b*){b_range}" if self.apply_strand_first
                  else f"(b{b_range})*")
        length = f"{'-' if self.length is None else self.length} bp"
        diff = f"{'-' if self.diff is None else round(100 * self.diff, 2)} %diff"
        return f"{a_repr} - {b_repr} ({length}, {diff})"

    @property
    def b_start_forward(self) -> Optional[int]:
        return (self.b_start if self.strand == 0 or not self.apply_strand_first
                else self.b_len - self.b_end if (self.b_len is not None
                                                 and self.b_end is not None)
                else None)

    @property
    def b_end_forward(self) -> Optional[int]:
        return (self.b_end if self.strand == 0 or not self.apply_strand_first
                else self.b_len - self.b_start if (self.b_len is not None
                                                   and self.b_start is not None)
                else None)

    def show(self, a_seq: Optional[str] = None, b_seq: Optional[str] = None,
             width: int = 100, twine_plot: bool = False):
        """Print the pairwise alignment like BLAST or as a "twine plot".
        `self.[a|b]_seq` are regarded as DNA sequences if `[a|b]_seq` are None.
        """
        assert self.cigar is not None, "`cigar` must be stored"
        if a_seq is None:
            assert isinstance(self.a_seq, str), "`self.a_seq` must be str"
            a_seq = self.a_seq
        if b_seq is None:
            assert isinstance(self.b_seq, str), "`self.b_seq` must be str"
            b_seq = self.b_seq

        # Cut out the subsequences of the alignment region
        a_seq = a_seq[self.a_start:self.a_end]
        b_seq = (b_seq[self.b_start:self.b_end] if self.strand == 0
                 else revcomp_seq(b_seq)[self.b_start:self.b_end] if self.apply_strand_first
                 else revcomp_seq(b_seq[self.b_start:self.b_end]))

        # Insert gaps to the subsequences
        fcigar = self.cigar.flatten()
        a_str, b_str = '', ''
        a_pos, b_pos = 0, 0
        for c in fcigar:
            if c == '=' or c == 'X':
                a_str += a_seq[a_pos]
                b_str += b_seq[b_pos]
                a_pos += 1
                b_pos += 1
            elif c == 'D':
                a_str += '-'
                b_str += b_seq[b_pos]
                b_pos += 1
            else:
                a_str += a_seq[a_pos]
                b_str += '-'
                a_pos += 1
        assert a_pos == len(a_seq) and b_pos == len(b_seq), \
            "Invalid CIGAR string for the alignment"

        if twine_plot:
            a_str, b_str, fcigar = map(lambda x: ''.join(x),
                                       zip(*[(' ', ' ', a_str[i]) if c == '='
                                             else (a_str[i], b_str[i], ' ')
                                             for i, c in enumerate(fcigar)]))
        a_str, b_str, fcigar = map(lambda x: split_seq(x, width),
                                   (a_str, b_str, fcigar))

        n_digit = max(map(lambda x: len(str(x)),
                          (self.a_start or 0,
                           self.a_end or 0,
                           self.b_start_forward or 0,
                           self.b_end_forward or 0)))
        a_pos = self.a_start
        b_pos = self.b_start_forward if self.strand == 0 else self.b_end_forward - 1

        for i in range(len(fcigar)):
            print(f"a:{'-' if a_pos is None else a_pos:{n_digit}}  {a_str[i]}")
            print(f"  {' ' * n_digit}  {fcigar[i]}")
            print(f"b:{'-' if b_pos is None else b_pos:{n_digit}}  {b_str[i]}")
            print("")
            if a_pos is not None:
                a_pos += len(a_str[i].replace('-', ''))
            if b_pos is not None:
                if self.strand == 0:
                    b_pos += len(b_str[i].replace('-', ''))
                else:
                    b_pos -= len(b_str[i].replace('-', ''))


EDLIB_MODE = {"global": "NW", "glocal": "HW", "prefix": "SHW"}


@dataclass(repr=False, eq=False)
class EdlibRunner:
    """Utility for running edlib with two sequences. Only "path" task is supported.

    usage:
      > er_global = EdlibRunner("global")
      > aln = er_global.align("acgtac", "accgac")
      > aln
      q[None:None] vs t[0:6]   (7 bp, 28.57% diff)

    positional arguments:
      @ mode : Must be one of {"global", "glocal", "prefix"}.

    optional arguments:
      @ revcomp      : If True, find reverse complement alignment as well.
      @ cyclic       : If True, perform cyclic alignment (`mode` must be `global`).
    """
    mode: str   # TODO: "proper" mode?
    revcomp: bool = True
    cyclic: bool = False

    def __post_init__(self):
        assert self.mode in EDLIB_MODE, f"Invalid mode: {self.mode}"
        assert not self.cyclic or self.mode == "global", \
            "Only global mode is supported for cyclic alignment"

    def align(self, query: str, target: str) -> Alignment:
        """Align `query` to `target`."""
        return (self._align_cyclic if self.cyclic
                else self._align)(query, target)

    def _align(self, query: str, target: str) -> Alignment:
        """Find best alignment with diff, cosidering strand if needed."""
        aln = self._run_edlib(query, target, 0)
        if self.revcomp:
            aln2 = self._run_edlib(query, target, 1)
            if aln.diff > aln2.diff:
                aln = aln2
        return aln

    def _run_edlib(self, query: str, target: str, strand: int) -> Alignment:
        """Run edlib with the two sequences.
        Multiple alignments with the same edit distance can exist, but here report only
        the first one.
        """
        e = edlib.align(query, target if strand == 0 else revcomp_seq(target),
                        mode=EDLIB_MODE[self.mode], task="path")
        start, end = e["locations"][0]
        cigar = Cigar(e["cigar"])
        return Alignment(a_seq=query, b_seq=target, strand=strand,
                         a_start=0, a_end=len(query), a_len=len(query),
                         b_start=start, b_end=end + 1, b_len=len(target),
                         length=cigar.aln_length,
                         diff=e["editDistance"] / cigar.aln_length,
                         cigar=cigar)

    def _align_cyclic(self, query: str, target: str) -> Alignment:
        """Map `query` to a duplicated `target` as a surrogate of cyclic alignment."""
        self.mode = "glocal"
        aln = self._align(query, target * 2)
        self.mode = "global"

        if aln.strand == 1:
            target = revcomp_seq(target)

        if aln.t_end > len(target):   # adjust the positions to the original coordinates
            aln.t_end -= len(target)
            if aln.t_start >= len(target):
                aln.t_start -= len(target)
        if aln.t_end != aln.t_start:   # adjust start/end positions if the mapping is not global
            # If the mapping is short, set `aln.t_start` to start=end position and re-calculate alignment.
            aln_s = self._realign_cyclic(query, target, aln.t_start)
            # If the mapping is redundant, check which of `aln.t_start` and `aln.t_end` is better for
            # start=end position of the alignment
            if aln.t_start < aln.t_end:
                aln_e = self._realign_cyclic(query, target, aln.t_end)
                if aln_s.diff > aln_e.diff:
                    aln_s = aln_e
            aln = aln_s

        return aln

    def _realign_cyclic(self, query: str, target: str, start_pos: int) -> Alignment:
        """Realign `query` globally assuming `target` starts from and end at `start_pos` cyclically."""
        aln = self._run_edlib(
            query, target[start_pos:] + target[:start_pos], 0)
        aln.t_start = start_pos if start_pos != len(target) else 0
        aln.t_end = start_pos if start_pos != 0 else len(target)
        return aln
