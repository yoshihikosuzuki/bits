from typing import NamedTuple, List, Dict
from collections import Counter, defaultdict
from BITS.seq.align import EdlibRunner

er_global = EdlibRunner("global", revcomp=False)


class DiscrepantSite(NamedTuple):
    seed_pos: int
    extra_pos: int
    edit_op: str
    base: str


def count_discrepant_sites(seed: str,
                           seqs: List[str]) -> Dict[DiscrepantSite, int]:
    """Count all discrepant sites between `seed` and each of `seqs`."""
    assert all([len(s) > 0 for s in seqs + [seed]]), "No empty strings"
    # TODO: how to decide "same variant?" especially for multiple variations
    # on same position (but slightly different among units)?
    discrepants = Counter()
    for seq in seqs:
        aln = er_global(seq, seed)
        gapped_seq, _ = aln.gapped_aligned_seq()
        seed_pos = 0
        extra_pos = 0   # positive values for continuous insertions
        for i, c in enumerate(aln.fcigar):   # fcigar: seed -> seq
            if c == '=':
                extra_pos = 0
            elif c == 'I':
                extra_pos += 1
            if c != '=':
                # TODO: multiple D on the same pos are aggregated
                discrepants[DiscrepantSite(seed_pos=seed_pos,
                                           extra_pos=extra_pos,
                                           edit_op=c,
                                           base=gapped_seq[i])] += 1
            if c != 'I':
                seed_pos += 1
        assert seed_pos == len(seed)
    return discrepants


def consensus_alt(in_seqs, seed_choice="original"):
    """Compute a consensus sequence among `seqs: List[str]` by a simple majority vote for each position
    of the alignment pileup that is made by globally aligning a seed sequence and each of the other sequences.
    """

    # Choose the seed and move it to the first element of `in_seqs`
    assert seed_choice in ("original", "median",
                           "longest"), "Invalid `seed_choise`"
    if seed_choice != "original":
        in_seqs = sorted(in_seqs, key=lambda x: len(x),
                         reverse=True)   # "longest"
        if seed_choice == "median":
            index = len(in_seqs) // 2
            in_seqs = [in_seqs[index]] + in_seqs[:index] + in_seqs[index + 1:]

    var_counts = count_discrepant_sites(in_seqs[0], in_seqs[1:])

    freqs = defaultdict(list)
    for (pos, subpos, _type, base), count in var_counts.items():
        freqs[(pos, subpos)].append((base, count))

    cons = ""

    for pos in range(len(in_seqs[0]) + 1):
        subpos = 1
        # insertions
        while (pos, subpos) in freqs:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in freqs[(pos, subpos)]:
                base_counts[base] = count
            base_counts['-'] = len(in_seqs) - sum(base_counts.values())
            cons += base_counts.most_common()[0][0]
            subpos += 1
        if pos == len(in_seqs[0]):
            break
        # others
        if (pos, 0) not in freqs:
            cons += in_seqs[0][pos]
        else:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in freqs[(pos, 0)]:
                base_counts[base] = count
            base_counts[in_seqs[0][pos]] = len(
                in_seqs) - sum(base_counts.values())
            cons += base_counts.most_common()[0][0]

    return cons.replace('-', '')
