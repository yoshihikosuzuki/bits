from typing import List, Optional, Sequence, Tuple, Type, Union

import pysam
from logzero import logger
from pyfastx import Fasta, Fastq

from ..util._proc import run_command
from ._type import BedRecord, FastaRecord, FastqRecord, SatRecord
from ._util import parse_region, split_seq


def _change_case(seq: str, case: str) -> str:
    assert case in (
        "original",
        "lower",
        "upper",
    ), "`case` must be 'original', 'lower', or 'upper'"
    return (
        seq if case == "original" else seq.lower() if case == "lower" else seq.upper()
    )


def load_fastx(
    in_fname: str,
    id_range: Optional[Union[int, Tuple[int, int]]] = None,
    case: str = "original",
    verbose: bool = True,
) -> List[Union[FastaRecord, FastqRecord]]:
    """Utility function in case one doesn't know sequence type."""
    FASTA_SUFFIXES = [".fa", ".fna", ".fasta"]
    FASTQ_SUFFIXES = [".fq", ".fastq"]
    if in_fname.endswith(
        tuple(FASTA_SUFFIXES + [f"{suf}.gz" for suf in FASTA_SUFFIXES])
    ):
        return load_fasta(in_fname, id_range, case, verbose)
    elif in_fname.endswith(
        tuple(FASTQ_SUFFIXES + [f"{suf}.gz" for suf in FASTQ_SUFFIXES])
    ):
        return load_fastq(in_fname, id_range, case, verbose)
    else:
        assert False, f"Cannot guess file type: {in_fname}"


def load_fasta(
    in_fname: str,
    id_range: Optional[Union[int, Tuple[int, int]]] = None,
    case: str = "original",
    verbose: bool = True,
) -> List[FastaRecord]:
    """Load (specified range of) a fasta file. Gzipped files are OK.

    positional arguments:
      @ in_fname : Input fasta file name.

    optional arguments:
      @ id_range : 1-indexed read ID or tuple of read IDs to be read.
      @ case     : Of the sequence to be stored.
                   Must be one of {"original", "lower", "upper"}.
    """
    is_single = isinstance(id_range, int)
    if id_range is None:
        seqs = [
            FastaRecord(name=name, seq=_change_case(seq, case))
            for name, seq in Fasta(in_fname, build_index=False, full_name=True)
        ]
    else:
        if is_single:
            id_range = (id_range, id_range)
        command = f"seqkit range -w0 -r{':'.join(map(str, id_range))} {in_fname}"
        out = run_command(command).strip().split("\n")
        assert len(out) % 2 == 0
        seqs = [
            FastaRecord(name=out[i * 2][1:], seq=_change_case(out[i * 2 + 1], case))
            for i in range(len(out) // 2)
        ]
    if verbose:
        logger.info(f"{in_fname}: {len(seqs)} sequences loaded")
    return seqs if not is_single else seqs[0]


def load_fastq(
    in_fname: str,
    id_range: Optional[Union[int, Tuple[int, int]]] = None,
    case: str = "original",
    verbose: bool = True,
) -> List[FastqRecord]:
    """Load (specified range of) a fastq file. Gzipped files are OK.

    positional arguments:
      @ in_fname : Input fastq file name.

    optional arguments:
      @ id_range : 1-indexed read ID or tuple of read IDs to be read.
      @ case     : Of the sequence to be stored.
                   Must be one of {"original", "lower", "upper"}.
    """
    is_single = isinstance(id_range, int)
    if id_range is None:
        seqs = [
            FastqRecord(name=name, seq=_change_case(seq, case), qual=qual)
            for name, seq, qual in Fastq(in_fname, build_index=False, full_name=True)
        ]
    else:
        if is_single:
            id_range = (id_range, id_range)
        command = f"seqkit range -r{':'.join(map(str, id_range))} {in_fname}"
        out = run_command(command).strip().split("\n")
        assert len(out) % 4 == 0
        seqs = [
            FastqRecord(
                name=out[i * 4][1:],
                seq=_change_case(out[i * 4 + 1], case),
                qual=out[i * 4 + 3],
            )
            for i in range(len(out) // 4)
        ]
    if verbose:
        logger.info(f"{in_fname}: {len(seqs)} sequences loaded")
    return seqs if not is_single else seqs[0]


def save_fasta(
    seqs: Union[FastaRecord, Sequence[FastaRecord]],
    out_fname: str,
    width: int = -1,
    verbose: bool = True,
) -> None:
    """If `width` > 0, newlines are inserted at every `width` bp."""
    assert width != 0, "`width` must not be 0"
    with open(out_fname, "w") as f:
        for seq in [seqs] if isinstance(seqs, FastaRecord) else seqs:
            _seq = seq.seq if width < 0 else "\n".join(split_seq(seq.seq, width))
            f.write(f">{seq.name}\n{_seq}\n")
    n_seq = len(seqs) if hasattr(seqs, "__len__") else 1
    if verbose:
        logger.info(f"{out_fname}: {n_seq} sequences saved")


def save_fastq(
    seqs: Union[FastqRecord, Sequence[FastqRecord]],
    out_fname: str,
    verbose: bool = True,
) -> None:
    with open(out_fname, "w") as f:
        for seq in [seqs] if isinstance(seqs, FastqRecord) else seqs:
            f.write(f"@{seq.name}\n{seq.seq}\n+\n{seq.qual}\n")
    n_seq = len(seqs) if hasattr(seqs, "__len__") else 1
    if verbose:
        logger.info(f"{out_fname}: {n_seq} sequences saved")


def load_bed(
    in_fname: str,
    attrs: List[Tuple[str, Type]] = [],
    verbose: bool = True,
) -> List[BedRecord]:
    """Load a bed file.

    Parameters
    ----------
    in_fname
        Input file name of a bed file
    attrs, optional
        A list of tuples of (attr_name, attr_type) for attributes after 4th column, by default []
    verbose, optional
        verbose mode, by default True

    Returns
    -------
        A list of records in the bed file
    """
    records = []
    with open(in_fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            data = line.strip().split("\t")
            if len(data) < 3:
                if verbose:
                    logger.warn(f"Ignoring record: {line.strip()}")
                continue
            r = BedRecord(chr=data[0], b=int(data[1]), e=int(data[2]))
            for (attr_name, attr_type), value in zip(attrs, data[3:]):
                r.__setattr__(attr_name, attr_type(value))
            records.append(r)
    if verbose:
        logger.info(f"{in_fname}: {len(records)} records loaded")
    return records


def load_trf(in_trf: str, verbose: bool = True) -> List[SatRecord]:
    sats = []
    with open(in_trf, "r") as f:
        for line in f:
            if line.startswith("@"):
                c = line.strip()[1:]
                continue
            data = line.strip().split()
            b = int(data[0])
            e = int(data[1])
            n_copy = float(data[3])
            unit_seq = data[13]
            sat = SatRecord(c, b, e, unit_seq, n_copy)
            sats.append(sat)
    if verbose:
        logger.info(f"{in_trf}: {len(sats)} records loaded")
    return sats


def load_bam(in_bam: str, region: Union[str, BedRecord]) -> List[pysam.AlignedSegment]:
    r = parse_region(region) if isinstance(region, str) else region
    return list(pysam.AlignmentFile(in_bam, "rb").fetch(r.chr, r.b, r.e))
