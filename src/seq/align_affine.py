import parasail
from BITS.seq.cigar import Cigar


def align_affine(query: str,
                 target: str,
                 match_score: int,
                 mismatch_cost: int,
                 gap_open_cost: int,
                 gap_extension_cost: int) -> Cigar:
    assert min(match_score,
               mismatch_cost,
               gap_open_cost,
               gap_extension_cost) >= 0, "Specify positive integers"
    return Cigar(parasail.nw_trace(query,
                                   target,
                                   gap_open_cost,
                                   gap_extension_cost,
                                   parasail.matrix_create("ACGT",
                                                          match_score,
                                                          -mismatch_cost))
                 .cigar.decode.decode('utf-8'))
