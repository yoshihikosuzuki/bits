from dataclasses import dataclass
from typing import List
from logzero import logger
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from .util import split_seq


@dataclass(frozen=True)
class FastaRecord:
    name: str
    seq: str


@dataclass(frozen=True)
class FastqRecord(FastaRecord):
    qual: str


def load_fasta(in_fname: str, case: str = "original") -> List[FastaRecord]:
    """`case` must be one of {"original" (default), "lower", "upper"}."""
    assert case in ("original", "lower", "upper"), \
        "`case` must be 'original', 'lower', or 'upper'"
    with open(in_fname, 'r') as f:
        seqs = [FastaRecord(name=name,
                            seq=(seq if case == "original"
                                 else seq.lower() if case == "lower"
                                 else seq.upper()))
                for name, seq in list(SimpleFastaParser(f))]
    logger.info(f"{len(seqs)} sequences loaded")
    return seqs


def load_fastq(in_fname: str, case: str = "original") -> List[FastqRecord]:
    """`case` must be one of {"original" (default), "lower", "upper"}."""
    assert case in ("original", "lower", "upper"), \
        "`case` must be 'original', 'lower', or 'upper'"
    with open(in_fname, 'r') as f:
        seqs = [FastqRecord(name=name,
                            seq=(seq if case == "original"
                                 else seq.lower() if case == "lower"
                                 else seq.upper()),
                            qual=qual)
                for name, seq, qual in FastqGeneralIterator(f)]
    logger.info(f"{len(seqs)} sequences loaded")
    return seqs


def save_fasta(reads: List[FastaRecord], out_fname: str, width: int = -1):
    """If `width` > 0, newlines are inserted at every `width` bp."""
    assert width != 0, "`width` must not be 0"
    with open(out_fname, 'w') as f:
        for read in reads:
            seq = read.seq if width < 0 else '\n'.join(split_seq(read.seq, width))
            f.write(f">{read.name}\n{seq}\n")


def save_fastq(reads: List[FastqRecord], out_fname: str):
    with open(out_fname, 'w') as f:
        for read in reads:
            f.write(f"@{read.name}\n{read.seq}\n+\n{read.qual}\n")
