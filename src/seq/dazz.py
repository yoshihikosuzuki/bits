from dataclasses import dataclass
from typing import Optional, Tuple, List
from logzero import logger
from BITS.util.proc import run_command
from .fastx import FastaRecord


@dataclass(eq=False)
class DazzRecord(FastaRecord):
    id: int

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"(id={repr(self.id)}, "
                f"name={repr(self.name)}, "
                f"seq=={repr(self.seq)})")


def load_db(db_fname: str,
            dbid_range: Optional[Tuple[int, int]] = None) -> List[DazzRecord]:
    if dbid_range is None:
        read_range = ""
        start_dbid = 1
    else:
        read_range = '-'.join(map(str, dbid_range))
        start_dbid = dbid_range[0]

    command = (f"DBshow {db_fname} {read_range} | "
               f"awk 'BEGIN {{first = 1}} "
               f"{{if (substr($0, 1, 1) == \">\") "
               f"{{if (first == 1) {{first = 0}} "
               f"else {{printf(\"%s\\t%s\\n\", header, seq)}}; "
               f"header = substr($0, 2); seq = \"\";}} "
               f"else {{seq = seq $0}}}} "
               f"END {{printf(\"%s\\t%s\\n\", header, seq)}}'")

    seqs = [DazzRecord(id=int(dbid), name=name, seq=seq)
            for dbid, (name, seq)
            in enumerate(map(lambda x: x.split('\t'),
                             run_command(command).strip().split('\n')),
                         start=start_dbid)]
    logger.info(f"{len(seqs)} sequences loaded")
    return seqs


def db_to_n_blocks(db_fname: str) -> int:
    """Extract the number of blocks from a DAZZ_DB file."""
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split('=')[1].strip())
    logger.error(f"No information on the number of blocks in {db_fname}")


def db_to_n_reads(db_fname: str) -> int:
    """Return the number of reads in a DAZZ_DB file."""
    return int(run_command(f"DBdump {db_fname} | awk 'NR == 1 {{print $3}}'").strip())
