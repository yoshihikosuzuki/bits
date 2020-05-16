import os
from dataclasses import dataclass, field
from typing import List
import numpy as np
from BITS.util.proc import run_command


@dataclass(repr=False, eq=False)
class ConsedRunner:
    """Class for running Consed especially with the options -V (variant call) and -G (variant graph), 
    which are not supported by the original consed package.
    Usage:
        r = ConsedRunner(variant_vector=True, variant_graph=True)
        r.run(["acgt", "acga", "acgat"])

    If you want to take consensus iteratively, do as follows:
        1) Take consensus iteratively using consed package, and
        2) Using the consensus sequence, call <run_consed>.
    """
    out_dir        : str   = "tmp"
    variant_vector : bool  = False
    variant_graph  : bool  = False
    min_var_frac   : float = 0.0   # minimum fraction of the variants to be shown
    display_width  : int   = 1000000
    seqs_fname     : str   = field(init=False)   # each instance has its own temp files
    out_fname      : str   = field(init=False) 

    def __post_init__(self):
        run_command(f"mkdir -p {self.out_dir}")
        self.seqs_fname = os.path.join(self.out_dir, f"in.seqs.{os.getpid()}")
        self.out_fname = os.path.join(self.out_dir, f"out.consed.{os.getpid()}")
        vv = "-V" if self.variant_vector else ""
        vg = f"-G{self.out_fname}" if self.variant_graph else ""
        vf = f"-t{self.min_var_frac}" if self.variant_vector or self.variant_graph else ""
        self.command = (f"consed {vv} {vg} {vf} -w{self.display_width} {self.seqs_fname} > {self.out_fname}")

    def _load_consensus(self):
        self.cons_seq = run_command(f"awk 'BEGIN {{seq = \"\"}} "
                                    f"$0 == \"\" {{exit}} {{seq = seq $0}} "
                                    f"END {{print seq}}' "
                                    f"{self.out_fname}").strip()

    def _load_variant_matrix(self):
        ret = run_command(f"grep -v '*' {self.out_fname} | "
                          f"grep -v 'versus' | "
                          f"sed -e '1,5d' | "
                          f"awk 'BEGIN {{i = 1}} "
                          f"{{if ($0 == \"\") {{i = 1}} else {{seq[i] = seq[i] $NF; i++}}}} "
                          f"END {{for (i = 1; i <= length(seq); i++) print seq[i]}}'").strip().split("\n")
        N = len(ret[0])   # number of input sequences
        M = len(ret)   # number of variants
        vmatrix = np.zeros((M, N), dtype=int)
        for i, r in enumerate(ret):
            vmatrix[i, :] = list(map(int, list(r.strip())))
        self.variant_matrix = vmatrix.T[1:-1, :]   # change to form of (sequences, variants)
        # TODO: implement extraction of variants information (move from src.old/encode.py)

    def run(self, in_seqs: List[str]):
        """Consensus sequence and variant matrix will be set to <self.cons_seq> and <self.variant_matrix>."""
        # Output <in_seqs> to a file
        with open(self.seqs_fname, "w") as f:
            for seq in in_seqs:
                f.write(f"{seq}\n")
        # Execute Consed
        run_command(self.command)
        # Load the results
        self._load_consensus()
        if self.variant_vector:
            self._load_variant_matrix()
