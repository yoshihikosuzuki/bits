import os
import re
from logzero import logger
import edlib

from .utils import run_command


edlib_mode = {"global": "NW", "glocal": "HW", "local": "SHW"}


def run_edlib(query, target, mode="global", swap_cigar=False):
    """
    Perform 1 vs 1 pairwise sequence alignment using edlib.
    You can specify one of the following modes: "global", "glocal", or "local".
    In the "glocal" mode, boundary-fixed <query> is mapped to <target>.

    Original CIGAR string output satisfies: "<query> == cigar(<target>[<start>:<end>])".
    You can swap the roles of <query> and <target> on the CIGAR using <swap_cigar=True> option (then "cigar(<query>) == <target>[<start>:<end>]").
    This is effective when you actually want to treat <query> as target (i.e. reference) in a SAM file output.
    """

    ret = edlib.align(query, target, mode=edlib_mode[mode], task="path")

    cigar = ret["cigar"]
    if swap_cigar is True:
        # Swap I and D in the cigar string in order to satisfy the definition of SAM format
        # when <query> is actually reference.
        # E.g. <query>: tandem repeat monomer, <target>: dimer
        cigar = cigar.replace("I", "?").replace("D", "I").replace("?", "D")

    # XXX: TODO: in "local" mode, we must eliminate in/dels in both sides (currently "diff" value in "local" mode is wrong)
    alignment_len = sum(list(map(int, re.findall(r'([0-9]+)', cigar))))   # including I/D
    diff = ret["editDistance"] / alignment_len   # must be 0 ~ 1
    start, end = ret["locations"][0]
    end += 1   # for consistency with that used in the fasta header

    return {"cigar": cigar,   # <query> == cigar(<target>[<start>:<end>]) if <swap_cigar> is False, otherwise cigar(<query>) == <target>[<start>:<end>]
            "diff": diff,   # percent difference within the overlapping interval
            "start": start,   # on <target>
            "end": end}


def run_consed(in_seqs, out_consensus, variant_vector=False, variant_graph=False, variant_fraction=0.3, display_width=80):
    """
    Take consensus of sequences using Consed.
    Variant vector and variant graph(s) will be output to <prefix of out_consensus>.V and .G, respectively.
    """
    
    logger.info("If Consed failed with 'Cannot align to first sequence', "
                "then do '$ sed -i -e \"<read_id + 1>d\" <aligned_units_fname>' and try again.")

    out_prefix = os.path.splitext(out_consensus)[0]

    vv = '-V' if variant_vector else ''
    vg = f"-G{out_prefix}" if variant_graph else ''
    vf = f"-t{variant_fraction}" if variant_vector or variant_graph else ''

    # Run Consed
    command = f"consed {vv} {vg} {vf} -w{display_width} {in_seqs}"
    consed_out = run_command(command)
    with open(f"{out_prefix}.consed", 'w') as out:
        out.write(f"{consed_out}")

    # Extract consensus sequence
    command = "awk 'BEGIN {seq = \"\"} $0 == \"\" {exit} {seq = seq $0} END {print seq}' %s" % f"{out_prefix}.consed"
    consensus_seq = run_command(command).strip()
    with open(out_consensus, 'w') as out:
        out.write(f">consensus/0/0_{len(consensus_seq)}\n{consensus_seq}\n")

    # Extract variant matrix
    command = "sed -e '0,/Variant/d' %s | sed -e '1,2d' | awk 'BEGIN {i = 1} {if ($0 == \"\") {i = 1} else {seq[i] = seq[i] $NF; i++}} END {for (i = 1; i <= length(seq); i++) print seq[i]}'" % f"{out_prefix}.consed"
    variant_matrix = run_command(command)
    with open(f"{out_prefix}.V", 'w') as out:
        out.write(f"{variant_matrix}")
