from typing import Union, Optional, Sequence, Tuple, List
from logzero import logger
from pyfastx import Fasta, Fastq
from ..util._proc import run_command
from ._type import FastaRecord, FastqRecord
from ._util import split_seq


def _change_case(seq: str, case: str) -> str:
    assert case in ("original", "lower", "upper"), \
        "`case` must be 'original', 'lower', or 'upper'"
    return (seq if case == "original"
            else seq.lower() if case == "lower"
            else seq.upper())


def load_fasta(in_fname: str,
               id_range: Optional[Union[int, Tuple[int, int]]] = None,
               case: str = "original",
               verbose: bool = True) -> List[FastaRecord]:
    """Load (specified range of) a fasta file.

    positional arguments:
      @ in_fname : Input fasta file name.

    optional arguments:
      @ id_range : 1-indexed read ID or tuple of read IDs to be read.
      @ case     : Of the sequence to be stored.
                   Must be one of {"original", "lower", "upper"}.
    """
    is_single = isinstance(id_range, int)
    if id_range is None:
        seqs = [FastaRecord(name=name,
                            seq=_change_case(seq, case))
                for name, seq in Fasta(in_fname, build_index=False, full_name=True)]
    else:
        if is_single:
            id_range = (id_range, id_range)
        command = f"seqkit range -w0 -r{':'.join(map(str, id_range))} {in_fname}"
        out = run_command(command).strip().split('\n')
        assert len(out) % 2 == 0
        seqs = [FastaRecord(name=out[i * 2][1:],
                            seq=_change_case(out[i * 2 + 1], case))
                for i in range(len(out) // 2)]
    if verbose:
        logger.info(f"{in_fname}: {len(seqs)} sequences loaded")
    return seqs if not is_single else seqs[0]


def load_fastq(in_fname: str,
               id_range: Optional[Union[int, Tuple[int, int]]] = None,
               case: str = "original",
               verbose: bool = True) -> List[FastqRecord]:
    """Load (specified range of) a fastq file.

    positional arguments:
      @ in_fname : Input fastq file name.

    optional arguments:
      @ id_range : 1-indexed read ID or tuple of read IDs to be read.
      @ case     : Of the sequence to be stored.
                   Must be one of {"original", "lower", "upper"}.
    """
    is_single = isinstance(id_range, int)
    if id_range is None:
        seqs = [FastqRecord(name=name,
                            seq=_change_case(seq, case),
                            qual=qual)
                for name, seq, qual in Fastq(in_fname, build_index=False, full_name=True)]
    else:
        if is_single:
            id_range = (id_range, id_range)
        command = f"seqkit range -r{':'.join(map(str, id_range))} {in_fname}"
        out = run_command(command).strip().split('\n')
        assert len(out) % 4 == 0
        seqs = [FastqRecord(name=out[i * 4][1:],
                            seq=_change_case(out[i * 4 + 1], case),
                            qual=out[i * 4 + 3])
                for i in range(len(out) // 4)]
    if verbose:
        logger.info(f"{in_fname}: {len(seqs)} sequences loaded")
    return seqs if not is_single else seqs[0]


def save_fasta(reads: Union[FastaRecord, Sequence[FastaRecord]],
               out_fname: str,
               width: int = -1,
               verbose: bool = True) -> None:
    """If `width` > 0, newlines are inserted at every `width` bp."""
    assert width != 0, "`width` must not be 0"
    with open(out_fname, 'w') as f:
        for read in ([reads] if isinstance(reads, FastaRecord) else reads):
            seq = (read.seq if width < 0
                   else '\n'.join(split_seq(read.seq, width)))
            f.write(f">{read.name}\n{seq}\n")
    n_seq = len(reads) if hasattr(reads, '__len__') else 1
    if verbose:
        logger.info(f"{out_fname}: {n_seq} sequences saved")


def save_fastq(reads: Union[FastqRecord, Sequence[FastqRecord]],
               out_fname: str,
               verbose: bool = True) -> None:
    with open(out_fname, 'w') as f:
        for read in ([reads] if isinstance(reads, FastqRecord) else reads):
            f.write(f"@{read.name}\n{read.seq}\n+\n{read.qual}\n")
    n_seq = len(reads) if hasattr(reads, '__len__') else 1
    if verbose:
        logger.info(f"{out_fname}: {n_seq} sequences saved")
