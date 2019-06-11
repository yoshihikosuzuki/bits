from logzero import logger
from BITS.util.proc import run_command


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
