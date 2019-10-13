from Bio import SeqIO
from .utils import split_seq


def load_fasta(in_fname):
    """Load a fasta file as a dictionary of {header: seq}.
    Dictionary is better than list beacuse you can use a header as a key.
    Dictionary is better than pd.Dataframe because it is simpler and can be easily converted into DF.
    """
    with open(in_fname, "r") as f:
        return dict(SeqIO.FastaIO.SimpleFastaParser(f))


def load_fastq(in_fname, only_qual=False):
    """Load a fastq file as a dictionary of {header: (seq, qual)},
    where <qual> is List[int], QV integer list.
    If <only_qual> is True, it instead returns {header: qual}.
    NOTE: this is slow due to Bio.SeqIO.parse, but I'm lazy.
    """
    if only_qual:
        return {record.name: record.letter_annotations["phred_quality"]
                for record in SeqIO.parse(in_fname, "fastq")}
    else:
       return {record.name: (record.seq, record.letter_annotations["phred_quality"])
                for record in SeqIO.parse(in_fname, "fastq")}


def save_fasta(reads, out_fname, sort=True, width=-1):
    """<reads> must be a dict of {header: seq}.
    If <sort> is True, the headers will be sorted.
    Newlines are inserted at every <width> bp in each sequence (-1 means no newlines).
    """
    with open(out_fname, "w") as f:
        for header, seq in sorted(reads.items()) if sort else reads.items():
            if width > 0:
                seq = "\n".join(split_seq(seq, width))
            f.write(f">{header}\n{seq}\n")
