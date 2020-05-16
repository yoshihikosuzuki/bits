from logzero import logger
from BITS.util.proc import run_command


def db_to_n_blocks(db_fname):
    """Extract the number of blocks from the db file."""
    with open(db_fname, 'r') as f:
        for line in f:
            if line.startswith("blocks"):
                return int(line.split('=')[1].strip())
    logger.error(f"No information on the number of blocks in {db_fname}")


def db_to_n_reads(db_fname):
    """Calculate the number of reads in `db_fname`."""
    return int(run_command(f"DBdump {db_fname} | awk 'NR == 1 {{print $3}}'").strip())
