from dataclasses import dataclass
from typing import Tuple
from logzero import logger
import edlib
from ._util import revcomp_seq, split_seq
from ._cigar import Cigar


@dataclass
class EdlibAlignment:
    """An alignment between two sequences computed by EdlibRunner.

    NOTE: `a_seq[a_start:a_end] == cigar(strand(b_seq[b_start:b_end]))`
          That is:
          1. `[a|b]_[start|end]` are positions on FORWARD sequence.
          2. `cigar` converts `b_aligned_seq` into `a_aligned_seq`.

    NOTE: Only `b_seq` can be cyclically aligned. `a_seq` never.
          Currently only global cyclic alignment is supported by EdlibRunner,
          meaning `b_start == b_end` indicates the alignment is cyclic.

    positional arguments:
      @ [a|b]_seq          : Input DNA sequenecs.
      @ strand             : 0 or 1.
      @ [a|b]_[start|end]  : Start/end position of the alignment.
      @ diff               : Sequence dissimilarity of the alignment.
      @ cigar              : Cigar string of the alignment.
    """
    a_seq: str
    b_seq: str
    strand: int
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    diff: float
    cigar: Cigar

    def __repr__(self) -> str:
        a_repr = f"a_seq[{self.a_start}:{self.a_end}]"
        b_repr = (f"b_seq[{self.b_start}:{self.b_end}]" if not self.is_cyclic
                  else f"b_seq[{self.b_start}:] + b_seq[:{self.b_end}]")
        if self.strand == 1:
            b_repr = f"*({b_repr})"
        return (f"{a_repr} ~ {b_repr} "
                f"({self.length} bp, {100 * self.diff:.2f} %diff)")

    @property
    def length(self) -> int:
        return self.cigar.aln_length

    @property
    def is_cyclic(self) -> bool:
        """Assuming only global alignment for cyclic alignment."""
        return self.b_start == self.b_end

    @property
    def a_aligned_seq(self) -> str:
        return self.a_seq[self.a_start:self.a_end]

    @property
    def b_aligned_seq(self) -> str:
        seq = (self.b_seq[self.b_start:self.b_end] if not self.is_cyclic
               else self.b_seq[self.b_start:] + self.b_seq[:self.b_end])
        return seq if self.strand == 0 else revcomp_seq(seq)

    def gapped_aligned_seqs(self) -> Tuple[str, str]:
        """Compute aligned sequences with explicit gap symbols inserted."""
        a_seq, b_seq = self.a_aligned_seq, self.b_aligned_seq
        a_gapped_seq, b_gapped_seq = '', ''
        a_pos, b_pos = 0, 0
        fcigar = self.cigar.flatten()
        for c in fcigar:
            if c != 'D':
                a_gapped_seq += a_seq[a_pos]
                a_pos += 1
            else:
                a_gapped_seq += '-'
            if c != 'I':
                b_gapped_seq += b_seq[b_pos]
                b_pos += 1
            else:
                b_gapped_seq += '-'
        assert a_pos == len(a_seq) and b_pos == len(b_seq), \
            "Invalid CIGAR string"
        assert len(a_gapped_seq) == len(b_gapped_seq) == len(fcigar), \
            "Inconsistent length"
        return a_gapped_seq, b_gapped_seq

    def show(self,
             a_name: str = "a",
             b_name: str = "b",
             width: int = 100,
             twist_plot: bool = False):
        """Print the pairwise alignment like BLAST or as a "twist plot"."""
        a_str, b_str = self.gapped_aligned_seqs()
        fcigar = self.cigar.flatten()
        if twist_plot:
            # Collapse bases on matches
            a_str, b_str, fcigar = map(lambda x: ''.join(x),
                                       zip(*[(' ', ' ', a_str[i]) if c == '='
                                             else (a_str[i], b_str[i], ' ')
                                             for i, c in enumerate(fcigar)]))
        a_str, b_str, fcigar = map(lambda x: split_seq(x, width),
                                   (a_str, b_str, fcigar))
        L_pos = max(map(lambda x: len(str(x)),
                        (self.a_start, self.a_end, self.b_start, self.b_end)))
        L_name = max(len(a_name), len(b_name))
        a_pos = self.a_start
        b_pos = self.b_start if self.strand == 0 else self.b_end - 1
        for i in range(len(fcigar)):
            print(f"{a_name:>{L_name}}:{a_pos:{L_pos}}  {a_str[i]}")
            print(f"{'':>{L_name}} {' ' * L_pos}  {fcigar[i]}")
            print(f"{b_name:>{L_name}}:{b_pos:{L_pos}}  {b_str[i]}")
            print("")
            a_pos += len(a_str[i].replace('-', ''))
            if self.strand == 0:
                b_pos += len(b_str[i].replace('-', ''))
            else:
                b_pos -= len(b_str[i].replace('-', ''))


EDLIB_MODE = {"global": "NW", "glocal": "HW", "prefix": "SHW"}


@dataclass(repr=False, eq=False)
class EdlibRunner:
    """Utility for running edlib with two sequences.

    NOTE: Only "path" task is currently supported.

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
    mode: str
    revcomp: bool = True
    cyclic: bool = False

    def __post_init__(self):
        assert self.mode in EDLIB_MODE, f"Invalid mode: {self.mode}"
        assert not self.cyclic or self.mode == "global", \
            "Currently only global mode is supported for cyclic alignment"
        assert not (self.revcomp and self.mode == "prefix"), \
            "Reverse complement cannot be cosidered for prefix alignment"

    def align(self,
              query: str,
              target: str) -> EdlibAlignment:
        """Align `query` to `target`."""
        if self.mode in ("glocal", "prefix") and len(query) > 1.1 * len(target):
            logger.warning(f"query ({len(query)} bp) is much longer than "
                           f"target ({len(target)} bp) unexpectedly")
        return (self._align_cyclic if self.cyclic
                else self._align)(query, target)

    def _align(self,
               query: str,
               target: str) -> EdlibAlignment:
        """Find best alignment with diff, cosidering strand if needed."""
        aln = self._run_edlib(query, target, strand=0)
        if self.revcomp:
            aln2 = self._run_edlib(query, target, strand=1)
            if aln.diff > aln2.diff:
                aln = aln2
        return aln

    def _run_edlib(self,
                   query: str,
                   target: str,
                   strand: int) -> EdlibAlignment:
        aln = edlib.align(query,
                          target if strand == 0 else revcomp_seq(target),
                          mode=EDLIB_MODE[self.mode],
                          task="path")
        # NOTE: multiple locations are possible, but here pick only the first one
        start, end = aln["locations"][0]
        end += 1
        if strand == 1:
            # Convert to positions on forward sequence
            start, end = len(target) - end, len(target) - start
        cigar = Cigar(aln["cigar"])
        return EdlibAlignment(a_seq=query,
                              b_seq=target,
                              strand=strand,
                              a_start=0,
                              a_end=len(query),
                              b_start=start,
                              b_end=end,
                              diff=aln["editDistance"] / cigar.aln_length,
                              cigar=cigar)

    def _align_cyclic(self,
                      query: str,
                      target: str) -> EdlibAlignment:
        """Map `query` to a duplicated `target` as a surrogate of cyclic alignment."""
        self.mode = "glocal"
        aln = self._align(query, target * 2)
        self.mode = "global"
        # Convert positions on duplicated sequence to those on original sequence
        aln.b_seq = target
        if aln.b_end > len(target):
            aln.b_end -= len(target)
            if aln.b_start >= len(target):
                aln.b_start -= len(target)
        if aln.b_end != aln.b_start:
            # Check which position is better for boundary of cyclic alignment
            aln_s = self._run_edlib_cyclic(query, target,
                                           aln.b_start, aln.strand)
            if aln.b_start < aln.b_end:
                aln_e = self._run_edlib_cyclic(query, target,
                                               aln.b_end, aln.strand)
                if aln_s.diff > aln_e.diff:
                    aln_s = aln_e
            aln = aln_s
        return aln

    def _run_edlib_cyclic(self,
                          query: str,
                          target: str,
                          t_boundary: int,
                          strand: int) -> EdlibAlignment:
        aln = self._run_edlib(query,
                              target[t_boundary:] + target[:t_boundary],
                              strand)
        aln.b_seq = target
        aln.b_start = aln.b_end = t_boundary
        return aln
