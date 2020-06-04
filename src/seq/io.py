from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, List, Dict
import numpy as np
from logzero import logger
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from BITS.util.type import ExplicitRepr
from BITS.util.proc import run_command
from .util import split_seq, ascii_to_phred


@dataclass
class FastaRecord(ExplicitRepr):
    """Sequence with name."""
    name: str
    seq: str

    def __repr__(self):
        return self._order_repr(["name", "seq"])

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass
class FastqRecord(FastaRecord):
    """Sequence with name and base qualities."""
    qual: str

    def __repr__(self):
        return self._order_repr(["name", "seq", "qual"])

    @property
    def qual_phred(self) -> np.ndarray:
        return np.array(list(map(ascii_to_phred, self.qual)),
                        dtype=np.int8)


@dataclass
class DazzRecord(FastaRecord):
    """Sequence with name and DAZZ_DB ID."""
    id: int

    def __repr__(self):
        return self._order_repr(["id", "name", "seq"])


@dataclass
class SeqInterval:
    """Class for an interval on a sequence."""
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


def load_fasta(in_fname: str,
               case: str = "original") -> List[FastaRecord]:
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


def load_fastq(in_fname: str,
               case: str = "original") -> List[FastqRecord]:
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


def save_fasta(reads: Sequence[FastaRecord],
               out_fname: str,
               width: int = -1):
    """If `width` > 0, newlines are inserted at every `width` bp."""
    assert width != 0, "`width` must not be 0"
    with open(out_fname, 'w') as f:
        for read in reads:
            seq = (read.seq if width < 0
                   else '\n'.join(split_seq(read.seq, width)))
            f.write(f">{read.name}\n{seq}\n")


def save_fastq(reads: Sequence[FastqRecord],
               out_fname: str):
    with open(out_fname, 'w') as f:
        for read in reads:
            f.write(f"@{read.name}\n{read.seq}\n+\n{read.qual}\n")


def load_db(db_fname: str,
            dbid_range: Optional[Tuple[int, int]] = None) -> List[DazzRecord]:
    """Load read IDs, original header names, and sequences from a DAZZ_DB file.
    `dbid_range` is e.g. `(1, 10)`, which is equal to `$ DBdump {db_fname} 1-10`.
    """
    n_reads = (db_to_n_reads(db_fname) if dbid_range is None
               else dbid_range[1] - dbid_range[0] + 1)
    seqs = [None] * n_reads
    i = 0
    command = (f"DBdump -rhs {db_fname} "
               f"{'' if dbid_range is None else '-'.join(map(str, dbid_range))}")
    for line in run_command(command).strip().split('\n'):
        line = line.strip()
        if line.startswith('R'):
            _, dazz_id = line.split()
        elif line.startswith('H'):
            _, _, prolog = line.split()
        elif line.startswith('L'):
            _, well, start, end = line.split()
            name = f"{prolog}/{well}/{start}_{end}"
        elif line.startswith('S'):
            _, _, seq = line.split()
            seqs[i] = DazzRecord(id=int(dazz_id), name=name, seq=seq)
            i += 1
    assert i == n_reads
    logger.info(f"{n_reads} sequences loaded")
    return seqs


def load_db_track(db_fname: str,
                  track_name: str,
                  dbid_range: Optional[Tuple[int, int]] = None) \
        -> Dict[int, List[SeqInterval]]:
    """Load track data of a DAZZ_DB file.
    Running a command like: `$ DBdump -r -m{track_name} {db_fname} {dbid_range}`.
    """
    tracks = {}
    count = 0
    command = (f"DBdump -r -m{track_name} {db_fname} "
               f"{'' if dbid_range is None else '-'.join(map(str, dbid_range))}")
    for line in run_command(command).strip().split('\n'):
        line = line.strip()
        if line.startswith('R'):
            _, read_id = line.split()
        elif line.startswith('T0'):
            poss = list(map(int, line.split()[2:]))
            assert len(poss) % 2 == 0, f"Cannot find pair of positions: {line}"
            tracks[int(read_id)] = [SeqInterval(start, end)
                                    for start, end in zip(poss[::2], poss[1::2])]
            count += len(poss) // 2
    logger.info(f"{count} intervals loaded from {len(tracks)} sequences.")
    return tracks


def db_to_n_blocks(db_fname: str) -> int:
    """Extract the number of blocks from a DAZZ_DB file."""
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split('=')[1].strip())
    logger.error(f"No information on the number of blocks in {db_fname}")


def db_to_n_reads(db_fname: str) -> int:
    """Return the number of reads in a DAZZ_DB file."""
    return int(run_command(f"DBdump {db_fname}").split('\n')[0].split()[-1])
