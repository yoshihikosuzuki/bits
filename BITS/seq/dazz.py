from BITS.util.proc import run_command


def db_to_n_reads(db_fname):
    """Calculate the number of reads in `db_fname`."""
    return int(run_command(f"DBdump {db_fname} | awk 'NR == 1 {{print $3}}'").strip())
