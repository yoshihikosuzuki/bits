from dataclasses import dataclass
from logzero import logger
import edlib
from .seq import revcomp
from .utils import run_command


@dataclass(repr=False)
class Cigar:
    string: str

    def iter(self):
        """
        Generator of a series of tuples (num. of bases, operator).
        """

        length = ""
        for c in self.string:
            if c == '=' or c == 'D' or c == 'I' or c == 'X' or c == 'N':
                yield (int(length), c)
                length = ""
            else:
                length += c

    @property
    def alignment_len(self):
        """
        Alignment length including gaps and masked regions.
        """

        return sum([l for l, op in self.iter()])

    def reversed(self):
        """
        Just reverse CIGAR without swapping in/del. Used for reverse complement.
        """

        return ''.join(reversed([f"{l}{op}" for l, op in self.iter()]))

    def flatten(self):
        """
        Convert CIGAR to a sequence of operations. Used for masking some specific intervals.
        """

        return FlattenCigar(''.join([op for l, op in self.iter() for i in range(l)]))

    def mask_intvl(self, intvl, ignore_op='D'):
        # TODO: how to do about I/D/X around intervals' boundaries

        cigar_f = self.flatten()
    
        starts = [i[0] for i in intvl]
        ends = [i[1] for i in intvl]
        index = 0
        pos = 0
        for i, c in enumerate(cigar_f):
            if index >= len(starts):
                break
            if i != 0 and c != ignore_op:
                pos += 1
            if pos > ends[index]:
                index += 1
            if index < len(starts) and starts[index] <= pos and pos <= ends[index]:   # NOTE: end-inclusive
                cigar_f[i] = 'N'

        return cigar_f.unflatten()


@dataclass(repr=False)
class FlattenCigar:
    string: str

    def unflatten(self):
        """
        Convert a series of oprations to a normal CIGAR string.
        """

        cigar = ""
        count = 0
        prev_c = self.string[0]
        for c in self.string:
            if c == prev_c:
                count += 1
            else:
                cigar += f"{count}{prev_c}"
                count = 1
                prev_c = c
        cigar += f"{count}{prev_c}"
        return Cigar(cigar)


edlib_mode = {"global": "NW", "glocal": "HW", "local": "SHW"}


class Align:
    def __init__(self, edlib_ret):
        # NOTE: coordinates are on the target sequence
        self.start, self.end = edlib_ret["locations"][0]
        self.end += 1   # for compatibility with python's slice

        self.edit_dist = edlib_ret["editDistance"]
        self.cigar = Cigar(edlib_ret["cigar"])
        self.length = self.cigar.alignment_len
        self.diff = self.edit_dist / self.length

        self.strand = 0   # set as forward (0) by default

    def show(self, show_cigar=False):
        logger.info(f"{'revcomp(' if self.strand == 1 else ''}target[{self.start}:{self.end}]{')' if self.strand == 1 else ''} ({self.length} bp, {self.diff:.3} diff)")
        if show_cigar:
            logger.info(f"{self.cigar.string}")


def _find_best_align(query,
                     target,
                     mode,
                     rc,
                     rc_pruning_diff_th,
                     strand_prior):
    if not rc:
        align = Align(edlib.align(query, target, mode=mode, task="path"))
    else:
        if strand_prior == 0:
            align_f = Align(edlib.align(query, target, mode=mode, task="path"))
            if align_f.diff < rc_pruning_diff_th:
                align = align_f
            else:
                target_rc = revcomp(target)
                align_rc = Align(edlib.align(query, target_rc, mode=mode, task="path"))
                if align_f.diff <= align_rc.diff:
                    align = align_f
                else:
                    align = align_rc
                    align.strand = 1
                    target = target_rc
        else:
            target_rc = revcomp(target)
            align_rc = Align(edlib.align(query, target_rc, mode=mode, task="path"))
            if align_rc.diff < rc_pruning_diff_th:
                align = align_rc
                align.strand = 1
                target = target_rc
            else:
                align_f = Align(edlib.align(query, target, mode=mode, task="path"))
                if align_f.diff <= align_rc.diff:
                    align = align_f
                else:
                    align = align_rc
                    align.strand = 1
                    target = target_rc

    return (align, target)


def run_edlib(query,
              target,
              mode,
              only_diff=False,
              rc=False,
              rc_pruning_diff_th=0.3,   # at most this diff if true strand
              strand_prior=0,   # 0 or 1
              cyclic=False,
              return_seq=False,   # in glocal alignment
              swap_cigar=False):
    """
    Perform 1 vs 1 pairwise sequence alignment using edlib.
    <mode> = {"global", "glocal", "local"}
    In "glocal" mode, boundary-fixed <query> is mapped to <target>.

    Originally, "<query> == cigar(<target>[<start>:<end>])" holds.
    If <swap_cigar> is True, "cigar(<query>) == <target>[<start>:<end>]" does.
    This is effective when you actually want to consider <query> as target
    (i.e. reference) in the SAM file format.

    NOTE: local alignment using edlib is not recommended.
    """

    assert mode in ["global", "glocal", "local"], f"Invalid mode name: {mode}"
    assert query != "" and target != "", "Input sequence must not be empty!"

    # NOTE: glocal cyclic possibly returns a shorter sequence than the original one
    global_cyclic = False
    if cyclic:
        target = target * 2
        if mode == "global":
            global_cyclic = True
            mode = "glocal"

    # Find best start/end and strand between the two sequences
    # If revcomp of target is the best, change target to it
    align, target = _find_best_align(query,
                                     target,
                                     edlib_mode[mode],
                                     rc,
                                     rc_pruning_diff_th,
                                     strand_prior)

    flag_cycle = False   # True if the alignment was actually cyclic
    if cyclic:
        # Adjust start/end points to the original coordinates
        half_pos = len(target) // 2
        if align.end > half_pos:
            align.end -= half_pos
            flag_cycle = not flag_cycle
            if align.start >= half_pos:
                align.start -= half_pos
                flag_cycle = not flag_cycle

        # Restore the duplicated target sequence to the original one
        target = target[:half_pos]

        # Adjust the difference between glo"b"al cyclic (expected) and glo"c"al mapping (firstly done)
        if global_cyclic:
            if flag_cycle:
                align.end = align.start
                target_re = target[align.start:] + target[:align.end]
            else:   # adjust so that the alignment is same as global mode
                align.start = 0
                align.end = half_pos
                target_re = target[align.start:align.end]
            # NOTE: taking alignment again is redundant, but very easy
            align_re = run_edlib(query, target_re, "global")
            align.diff = align_re.diff
            align.cigar = align_re.cigar

    if only_diff:
        return align.diff

    if return_seq:
        if flag_cycle:
            align.seq = target[align.start:] + target[:align.end]
        else:
            align.seq = target[align.start:align.end]

    if align.strand == 1:
        tmp = align.start
        align.start = len(target) - align.end
        align.end = len(target) - tmp

    if swap_cigar:
        align.cigar = align.cigar.replace("I", "?").replace("D", "I").replace("?", "D")

    return align


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
