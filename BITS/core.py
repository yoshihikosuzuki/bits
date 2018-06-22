import os
from logzero import logger
from .utils import run_command


def extract_adapters_from_bax(in_bax, out_adapters):
    """
    Extract all adapter sequences detected in PacBio reads
    """
    
    from pbcore.io import BasH5Reader

    with open(out_adapters, 'w') as out:
        with BasH5Reader(in_bax) as f:
            for r in f:
                for a in r.adapters:
                    out.write(a.basecalls() + '\n')


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


def consensus_adaptors(in_bax, out_adapters='adapters', out_consensus='adapter.consensus.fasta'):
    """
    Extract all adapter sequences from *.bax.h5, and then take consensus of them.
    """

    extract_adapters_from_bax(in_bax, out_adapters)
    run_consed(out_adapters, out_consensus)
