from dataclasses import dataclass
from logzero import logger
import edlib
from BITS.util.proc import run_command
from .core import revcomp, Cigar

edlib_mode = {"global": "NW", "glocal": "HW", "local": "SHW"}


@dataclass(repr=False, eq=False)
class Alignment:
    """
    Class for specifying the fields an alignment should have.
      query == <cigar>(<strand>(target[<t_start>:<t_end>]))
    """
    q_start : int   = None
    q_end   : int   = None
    t_start : int   = None
    t_end   : int   = None
    strand  : int   = 0      # 0 or 1
    length  : int   = None
    diff    : float = None
    cigar   : Cigar = None

    def __str__(self):
        return (f"q[{self.q_start}:{self.q_end}] vs "
                f"t{'' if self.strand == 0 else '*'}[{self.t_start}:{self.t_end}]   "
                f"({self.length} bp, {self.diff:.3} %diff)")

    def mapped_seq(self, target):
        """Returns the substring of <target> (as forward sequence) alignd to the query."""
        ret = (target[self.t_start:self.t_end] if self.t_start != self.t_end
               else target[self.t_start:] + target[self.t_start])   # accept only global for cyclic
        if self.strand == 1:
            ret = revcomp(ret)
        return ret


def _set_cigar(aln, e):
    if e["cigar"] is None:
        return
    aln.cigar = Cigar(e["cigar"])
    aln.length = aln.cigar.alignment_len
    aln.diff = e["editDistance"] / aln.length


def _realign_cyclic(start_pos, q, t, s):
    e = edlib.align(q, t[start_pos:] + t[:start_pos], "NW", task="path")
    a = Alignment(t_start=start_pos if start_pos != len(t) else 0,
                  t_end=start_pos if start_pos != 0 else len(t),
                  strand=s)
    _set_cigar(a, e)
    return a

        
@dataclass(repr=False, eq=False)
class EdlibRunner:
    """
    Class for running edlib for pairwise alignment of 2 sequences with several options.
    In "glocal" mode, <query> is mapped to <target>.
    Usage:
        r = EdlibRunner("global", rc=True, cyclic=False)
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
            if aln.t_end != aln.t_start:   # mapping doesn't or does redundantly cover <target>
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


def run_consed(in_fname,
               out_fname="out.consed",
               variant_vector=False,
               variant_graph=False,
               variant_fraction=0.3,
               display_width=1000000):
    """
    Only from file to file. No iterations.
    """

    vv = '-V' if variant_vector else ''
    vg = f"-G{out_fname}" if variant_graph else ''
    vf = f"-t{variant_fraction}" if variant_vector or variant_graph else ''
    run_command(f"consed {vv} {vg} {vf} -w{display_width} {in_fname} > {out_fname}")


def consed_to_consensus(in_consed):
    run_command(f"awk 'BEGIN {{seq = \"\"}} $0 == \"\" {{exit}} {{seq = seq $0}} END {{print seq}}' {in_consed} > {in_consed}.cons_seq")


def consed_to_varmat(in_consed):
    run_command(f"grep -v '*' {in_consed} | grep -v 'versus' | sed -e '1,5d' | awk 'BEGIN {{i = 1}} {{if ($0 == \"\") {{i = 1}} else {{seq[i] = seq[i] $NF; i++}}}} END {{for (i = 1; i <= length(seq); i++) print seq[i]}}' > {in_consed}.V")
