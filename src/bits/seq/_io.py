from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple, Type, Union

import pysam
from logzero import logger
from pyfastx import Fasta, Fastq

from ..util._proc import run_command
from ._type import BedRecord, FastaRecord, FastqRecord, GffRecord, SatRecord, SegRecord
from ._util import split_seq


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
    region: Optional[Union[str, SegRecord]] = None,
    attr_cols: Optional[List[int]] = None,
    attrs: List[Tuple[str, Type]] = [],
    by_chrom: bool = False,
    verbose: bool = True,
) -> Union[List[BedRecord], Dict[str, List[BedRecord]]]:
    """Load a bed file.

    Parameters
    ----------
    in_fname
        Input file name of a bed file
    attr_cols, optional
        A list of (1-based) indexes of columns for attributes, by default None
    attrs, optional
        A list of tuples of (attr_name, attr_type) for attributes after 4th column, by default []
    verbose, optional
        verbose mode, by default True

    Returns
    -------
        A list of records in the bed file
    """
    if region is not None and isinstance(region, str):
        region = SegRecord.from_string(region)
    if attr_cols is not None:
        assert len(attr_cols) == len(
            attrs
        ), f"Inconsistent len(attrs) {len(attrs)} vs len(attr_cols) {len(attr_cols)}"

    records = []
    with open(in_fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            data = line.strip().split("\t")
            if len(data) < 3:
                if verbose:
                    logger.warning(f"Ignoring record: {line.strip()}")
                continue
            r = BedRecord(chr=data[0], b=int(data[1]), e=int(data[2]))
            if region is not None and not (
                region.chr == r.chr
                and (region.b is None or region.b <= r.b)
                and (region.e is None or r.e <= region.e)
            ):
                continue
            for (attr_name, attr_type), value in zip(
                attrs,
                data[3:] if attr_cols is None else [data[n - 1] for n in attr_cols],
            ):
                r.__setattr__(attr_name, attr_type(value))
            records.append(r)
    if verbose:
        logger.info(f"{in_fname}: {len(records)} records loaded")
    if not by_chrom:
        return records
    else:
        records_by_chrom = defaultdict(list)
        for r in records:
            records_by_chrom[r.chr].append(r)
        return records_by_chrom


def filter_bed(
    data: Sequence[BedRecord],
    region: Optional[Union[str, SegRecord]],
    verbose: bool = True,
) -> List[BedRecord]:
    """Filter BedRecords that are already loaded."""
    if region is None:
        return data

    if isinstance(region, str):
        region = SegRecord.from_string(region)
    n_before = len(data)
    records = list(
        filter(
            lambda x: x.chr == region.chr and region.b <= x.b and x.e <= region.e,
            data,
        )
    )

    if verbose:
        logger.info(f"{n_before} -> {len(records)} records")
    return records


def load_trf(in_trf: str, verbose: bool = True) -> List[SatRecord]:
    sats = []
    with open(in_trf, "r") as f:
        for line in f:
            if line.startswith("@"):
                chrom = line.strip()[1:]
                continue
            data = line.strip().split()
            sat = SatRecord(
                chr=chrom,
                b=int(data[0]),
                e=int(data[1]),
                unit_seq=data[13],
                n_copy=float(data[3]),
            )
            sats.append(sat)
    if verbose:
        logger.info(f"{in_trf}: {len(sats)} records loaded")
    return sats


def load_bam(
    in_bam: str,
    region: Optional[Union[str, SegRecord]] = None,
    require_span: bool = False,
    verbose: bool = True,
) -> List[pysam.AlignedSegment]:
    if region is None:
        mappings = list(pysam.AlignmentFile(in_bam).fetch())
    else:
        mappings = list(
            pysam.AlignmentFile(in_bam).fetch(
                region=region.to_string() if isinstance(region, SegRecord) else region
            )
        )

    if require_span:
        assert region is not None, "`region` must be specified for `require_span`"
        if isinstance(region, str):
            region = SegRecord.from_string(region)
        assert (
            region.b is not None and region.e is not None
        ), "Both start/end positions must be specified in `region` for `require_span`"
        mappings = list(
            filter(
                lambda m: m.reference_start <= region.b - 1
                and m.reference_end >= region.e + 1,
                mappings,
            )
        )
    if verbose:
        logger.info(f"{in_bam}: {len(mappings)} records loaded")
    return mappings


def load_vcf(
    in_fname: str,
    region: Optional[Union[str, SegRecord]] = None,
    verbose: bool = True,
) -> List[pysam.VariantRecord]:
    """Load a vcf file.

    Parameters
    ----------
    in_fname
        Input file name of a vcf file
    verbose, optional
        verbose mode, by default True

    Returns
    -------
        A list of records in the vcf file
    """
    if region is None:
        variants = list(pysam.VariantFile(in_fname).fetch())
    else:
        variants = list(
            pysam.VariantFile(in_fname).fetch(
                region=region.to_string() if isinstance(region, SegRecord) else region
            )
        )

    if verbose:
        logger.info(f"{in_fname}: {len(variants)} records loaded")
    return variants


def load_gff(
    in_fname: str,
    region: Optional[Union[str, SegRecord]] = None,
    filter_type: Optional[str] = None,
    guess_attr_type: bool = True,
    attr_sep: Optional[str] = '=',
    verbose: bool = True,
) -> List[GffRecord]:
    """Load a gff file.

    Parameters
    ----------
    in_fname
        Input file name of a gff file
    region, optional
        Region to be loaded
    filter_type, optional
        Type of records to be loaded. e.g. "gene"
    guess_attr_type
        If True, automatically guess the type (only int and float) of each attribute
    """

    def _guess_type(v):
        try:
            return int(v)
        except ValueError:
            pass
        try:
            return float(v)
        except ValueError:
            pass
        return v

    records = []
    with open(in_fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            data = line.strip().split("\t")
            chrom, source, _type, b, e, _, strand, _ = data[:8]
            if filter_type is not None and _type != filter_type:
                continue
            b, e = int(b), int(e)
            if region is not None and not (
                region.chr == chrom and region.b <= b and e < region.e
            ):
                continue
            r = GffRecord(
                chr=chrom,
                b=b,
                e=e,
                forward=(strand == "+"),
                type=_type,
                source=source
            )
            if len(data) == 9:
                attrs = data[8]
                for k, v in map(lambda attr: attr.strip().split(attr_sep), attrs.strip(";").split(";")):
                    r.__setattr__(k, _guess_type(v) if guess_attr_type else v)
            elif len(data) >= 10:
                r.attrs = data[8:]
            records.append(r)
    if verbose:
        logger.info(f"{in_fname}: {len(records)} records loaded")
    return records
