from dataclasses import dataclass
from logzero import logger
import edlib
from .utils import revcomp
from .cigar import Cigar

edlib_mode = {"global": "NW", "glocal": "HW", "local": "SHW"}


@dataclass(repr=False, eq=False)
class Alignment:
    """query == <cigar>(<strand>(target[<t_start>:<t_end>]))"""
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
                f"({self.length} bp, {100 * self.diff:.3}% diff)")

    def mapped_seq(self, target):
        """Returns the substring of <target> (as forward sequence) alignd to the query."""
        ret = (target[self.t_start:self.t_end] if self.t_start != self.t_end
               else target[self.t_start:] + target[:self.t_start])   # accept only global for cyclic
        if self.strand == 1:
            ret = revcomp(ret)
        return ret


def _set_cigar(aln, edlib_ret):
    """Set CIGAR-related attributes in the given <aln> (assuming passed by reference)."""
    if edlib_ret["cigar"] is None:
        return
    aln.cigar = Cigar(edlib_ret["cigar"])
    aln.length = aln.cigar.alignment_len
    aln.diff = edlib_ret["editDistance"] / aln.length


def _realign_cyclic(start_pos, query, target, strand):
    """Realign query globally assuming <target> starts from and end at <start_pos> cyclically."""
    edlib_ret = edlib.align(query, target[start_pos:] + target[:start_pos], "NW", task="path")
    aln = Alignment(t_start=start_pos if start_pos != len(target) else 0,
                    t_end=start_pos if start_pos != 0 else len(target),
                    strand=strand)
    _set_cigar(aln, edlib_ret)
    return aln

        
@dataclass(repr=False, eq=False)
class EdlibRunner:
    """Class for running edlib for pairwise alignment of 2 sequences with several options.
    In "glocal" mode, <query> is mapped to <target>.
    Usage:
        r = EdlibRunner("global", revcomp=True, cyclic=False)
        a = r.align("acgtac", "accgac")
    """
    mode          : str
    revcomp       : bool  = True   # find reverse complement alignment as well if True
    cyclic        : bool  = False   # do cyclic alignment if True; currently only global cyclic is supported
    only_diff     : bool  = False   # do not trace back nor calculate CIGAR; a bit faster
    strand_prior  : int   = 0       # prior information on the strand [0 or 1]
    max_true_diff : float = 0.0     # skip alignment of the other strand if diff is already less than this

    def __post_init__(self):
        assert self.mode in edlib_mode.keys(), f"Invalid mode name: {self.mode}"
        if self.cyclic:
            assert self.mode == "global", "Only global mode is now supported for cyclic alignment"
            assert not self.only_diff, "Invalid combination of <cyclic> and <only_diff> options."
            self.mode = "glocal"
        if not self.revcomp and self.strand_prior == 1:
            logger.error("Invalid combination of <revcomp> and <strand_prior> options.")

    def _run(self, query, target, strand):
        e = edlib.align(query,
                        target if strand == 0 else revcomp(target),
                        mode=edlib_mode[self.mode],
                        task="distance" if self.only_diff else "path")
        aln = Alignment(t_start=e["locations"][0][0], t_end=e["locations"][0][1] + 1, strand=strand)
        _set_cigar(aln, e)
        return aln

    def _find_best_alignment(self, query, target):
        aln = self._run(query, target, self.strand_prior)
        if self.revcomp and aln.diff > self.max_true_diff:
            aln2 = self._run(query, target, 1 - self.strand_prior)
            if aln.diff > aln2.diff:
                aln = aln2
        return aln

    def align(self, query, target):
        """Take alignment between <query> and <target>."""
        if self.cyclic:
            target += target   # duplicate the target sequence
        aln = self._find_best_alignment(query, target)
        if aln.strand == 1:
            target = revcomp(target)
        if self.cyclic:
            target = target[:len(target) // 2]   # restore the target sequence
            if aln.t_end > len(target):   # adjust the positions to the original coordinates
                aln.t_end -= len(target)
                if aln.t_start >= len(target):
                    aln.t_start -= len(target)
            if aln.t_end != aln.t_start:   # mapping doesn"t or does redundantly cover <target>
                # If it is short mapping, then set it as <aln.t_start>.
                aln_s = _realign_cyclic(aln.t_start, query, target, aln.strand)
                # Choose better start(=end) position from <aln.t_start> and <aln.t_end> if redundant mapping.
                if aln.t_start < aln.t_end:
                    aln_e = _realign_cyclic(aln.t_end, query, target, aln.strand)
                    if aln_s.diff > aln_e.diff:
                        aln_s = aln_e
                aln = aln_s
        if aln.strand == 1:
            aln.t_start, aln.t_end = len(target) - aln.t_end, len(target) - aln.t_start
        return aln
