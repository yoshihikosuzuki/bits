from dataclasses import dataclass
import edlib
from .util import revcomp_seq
from .cigar import Cigar

edlib_mode = {"global": "NW", "glocal": "HW", "prefix": "SHW"}


@dataclass(repr=False, eq=False)
class Alignment:
    """query == `cigar`(`strand`(target[`t_start`:`t_end`]))"""
    q_start : int   = None
    q_end   : int   = None
    t_start : int   = None
    t_end   : int   = None
    strand  : int   = 0      # 0 or 1
    length  : int   = None
    diff    : float = None
    cigar   : Cigar = None

    def __repr__(self):
        return (f"q[{self.q_start}:{self.q_end}] vs "
                f"t{'' if self.strand == 0 else '*'}[{self.t_start}:{self.t_end}]   "
                f"({self.length} bp, {round(100 * self.diff, 2)}% diff)")

    def mapped_seq(self, target):
        """Returns the substring of `target` (as forward sequence) alignd to the query."""
        ret = (target[self.t_start:self.t_end] if self.t_start != self.t_end
               else target[self.t_start:] + target[:self.t_start])   # accept only global for cyclic
        if self.strand == 1:
            ret = revcomp_seq(ret)
        return ret


@dataclass(repr=False, eq=False)
class EdlibRunner:
    """Class for running edlib for pairwise alignment of 2 sequences with several options.
    The "distance" task of edlib is not supported because precise alignment diff requires precise
    alignment length, which is calculated from CIGAR not provided in the "distance" task.

    Usage example:
        > r = EdlibRunner("global", revcomp=True, cyclic=False)
        > a = r.align("acgtac", "accgac")

    positional arguments:
      @ mode <str> : Must be "global", "glocal" or "prefix".

    optional arguments:
      @ revcomp      <bool> [True]  : If `True`, find reverse complement alignment as well.
      @ cyclic       <bool> [True]  : If `True`, perform cyclic alignment. `mode` must be `global`.
      @ strand_prior <int>  [0]     : Prior information on the strand. Must be 0 or 1.
    """
    mode          : str
    revcomp       : bool  = True
    cyclic        : bool  = False
    strand_prior  : int   = 0

    def __post_init__(self):
        assert self.mode in edlib_mode.keys(), f"Invalid mode name: {self.mode}"
        assert self.revcomp or self.strand_prior == 0, "`strand_prior` must be 0 when `revcomp=False`"
        assert not self.cyclic or self.mode == "global", "Only global mode is allowed for cyclic alignment"

    def align(self, query, target):
        """Take alignment between `query` and `target`."""
        aln = (self._find_best_alignment(query, target) if not self.cyclic
               else self._align_cyclic(query, target))
        if aln.strand == 1:
            aln.t_start, aln.t_end = len(target) - aln.t_end, len(target) - aln.t_start
        return aln

    def _align_cyclic(self, query, target):
        """Map `query` to a duplicated `target` as a surrogate of cyclic alignment."""
        self.mode = "glocal"
        aln = self._find_best_alignment(query, target * 2)
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

    def _realign_cyclic(self, query, target, start_pos):
        """Realign `query` globally assuming `target` starts from and end at `start_pos` cyclically."""
        aln = self._run(query, target[start_pos:] + target[:start_pos], 0)
        aln.t_start = start_pos if start_pos != len(target) else 0
        aln.t_end = start_pos if start_pos != 0 else len(target)
        return aln

    def _find_best_alignment(self, query, target):
        aln = self._run(query, target, self.strand_prior)
        if self.revcomp:
            aln2 = self._run(query, target, 1 - self.strand_prior)
            if aln.diff > aln2.diff:
                aln = aln2
        return aln

    def _run(self, query, target, strand):
        # Run edlib
        e = edlib.align(query, target if strand == 0 else revcomp_seq(target),
                        mode=edlib_mode[self.mode], task="path")

        # Multiple alignments of the same edit distance can exist, but that is very low probability
        # and thus greedily choose the first one
        start, end = e["locations"][0]
        end += 1
        aln = Alignment(t_start=start, t_end=end, strand=strand)

        # Set CIGAR and alignment length
        aln.cigar = Cigar(e["cigar"])
        aln.length = aln.cigar.alignment_len
        aln.diff = e["editDistance"] / aln.length

        return aln
