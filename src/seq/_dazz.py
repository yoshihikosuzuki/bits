import os
from typing import Optional, Sequence, Tuple, List, Dict
from logzero import logger
from ..util import run_command
from ._seqtype import DazzRecord, SeqInterval


def fasta_to_db(fasta_fname: str,
                db_prefix: str,
                db_block_size: int = 500,
                db_type: str = "db") -> None:
    assert db_type in ("db", "dam"), \
        "`db_type` must be one of {'db', 'dam'}"
    run_command(f"fasta2{db_type.upper()} {db_prefix} {fasta_fname}")
    run_command(f"DBsplit -s{db_block_size} {db_prefix}")
    n_reads = db_to_n_reads(f"{db_prefix}.{db_type}")
    n_blocks = db_to_n_blocks(f"{db_prefix}.{db_type}")
    logger.info(f"{n_reads} reads and {n_blocks} blocks")


def load_db(db_fname: str,
            dbid_range: Optional[Tuple[int, int]] = None) -> List[DazzRecord]:
    """Load read IDs, original header names, and sequences from a DAZZ_DB file.
    `dbid_range` is e.g. `(1, 10)`, which is equal to `$ DBdump {db_fname} 1-10`.
    """
    mode = os.path.splitext(db_fname)[1]   # ".db" or ".dam"
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
            if mode == ".db":
                _, _, prolog = line.split()
            else:
                name = line.split('>')[1]
        elif line.startswith('L'):
            if mode == ".db":
                _, well, start, end = line.split()
                name = f"{prolog}/{well}/{start}_{end}"
        elif line.startswith('S'):
            _, _, seq = line.split()
            seqs[i] = DazzRecord(id=int(dazz_id), name=name, seq=seq)
            i += 1
    assert i == n_reads
    logger.info(f"{db_fname}: {n_reads} sequences loaded")
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
    logger.info(f"{db_fname} ({track_name}): "
                f"{count} intervals loaded from {len(tracks)} sequences.")
    return tracks


def db_to_n_blocks(db_fname: str) -> int:
    """Extract the number of blocks from a DAZZ_DB file."""
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split('=')[1].strip())
    logger.error(f"{db_fname}: No information on the number of blocks")


def db_to_n_reads(db_fname: str) -> int:
    """Return the number of reads in a DAZZ_DB file."""
    return int(run_command(f"DBdump {db_fname}").split('\n')[0].split()[-1])
