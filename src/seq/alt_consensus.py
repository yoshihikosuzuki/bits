from typing import NamedTuple, Sequence, List, Dict
from collections import Counter, defaultdict
from BITS.seq.align import EdlibRunner

er_global = EdlibRunner("global", revcomp=False)


class DiscrepantSite(NamedTuple):
    seed_pos: int
    extra_pos: int   # positions between adjacent seed positions
    edit_op: str
    base: str


def count_discrepants(seed: str,
                      seqs: Sequence[str]) -> Dict[DiscrepantSite, int]:
    """Count all discrepant sites between `seed` and each of `seqs`."""
    assert len(seed) > 0 and all([len(s) > 0 for s in seqs]), \
        "No empty strings"
    # TODO: how to decide "same variant?" especially for multiple variations
    # on same position (but slightly different among units)?
    discrepants = Counter()
    for seq in seqs:
        aln = er_global.align(seq, seed)
        gapped_seq, _ = aln.gapped_aligned_seqs()
        seed_pos = 0
        extra_pos = 0   # positive values for continuous insertions
        for i, c in enumerate(aln.cigar.flatten()):   # fcigar: seed -> seq
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


def consensus_alt(in_seqs: List[str],
                  seed_choice: str = "original") -> str:
    """Compute a consensus sequence by a simple majority vote for each position
    of the alignment pileup made from global alignments of sequqnces to a seed
    sequence."""
    assert seed_choice in ("original", "median", "longest"), \
        "`seed_choise` must be one of {'original', 'median', 'longest'}"
    if seed_choice != "original":
        in_seqs = sorted(in_seqs, key=lambda x: len(x), reverse=True)
        if seed_choice == "median":
            index = len(in_seqs) // 2
            in_seqs = [in_seqs[index]] + in_seqs[:index] + in_seqs[index + 1:]
    # Count discrepancies on the seed
    discrepant_counts = count_discrepants(in_seqs[0], in_seqs[1:])
    discrepants = defaultdict(list)
    for dis, count in discrepant_counts.items():
        discrepants[(dis.seed_pos, dis.extra_pos)].append((dis.base, count))
    # Majority vote for each position on the seed
    cons_seq = ""
    seed_seq = in_seqs[0]
    for pos in range(len(seed_seq) + 1):
        # Insertions
        extra_pos = 1
        while (pos, extra_pos) in discrepants:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in discrepants[(pos, extra_pos)]:
                base_counts[base] = count
            base_counts['-'] = len(in_seqs) - sum(base_counts.values())
            cons_seq += base_counts.most_common()[0][0]
            extra_pos += 1
        if pos == len(seed_seq):
            break
        # Others
        seed_base = seed_seq[pos]
        if (pos, 0) not in discrepants:
            cons_seq += seed_base
        else:
            base_counts = Counter({'a': 0, 'c': 0, 'g': 0, 't': 0, '-': 0})
            for base, count in discrepants[(pos, 0)]:
                base_counts[base] = count
            base_counts[seed_base] = len(in_seqs) - sum(base_counts.values())
            cons_seq += base_counts.most_common()[0][0]
    return cons_seq.replace('-', '')
