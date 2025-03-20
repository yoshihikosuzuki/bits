from dataclasses import InitVar, dataclass, field
from typing import Optional

import igv_notebook
from logzero import logger


@dataclass(repr=False, eq=False)
class IGVbrowser:
    """Wrapped IGV browser object.

    For details, see:
        - https://github.com/igvteam/igv-notebook
        - https://github.com/igvteam/igv.js/wiki/Tracks-2.0

    Usage example in Jupyter:
      > IGVbrowser("/path/to/fasta", region="chr1:10000-20000") \\
      >     .add_bam("/path/to/bam", track_name="bam")
    """

    ref_fname_or_genome_name: InitVar[str]
    ref_id: InitVar[Optional[str]] = None
    ref_name: InitVar[Optional[str]] = "ref"
    region: InitVar[Optional[str]] = None

    b: igv_notebook.Browser = field(init=False)

    def __post_init__(self, ref_fname_or_genome_name, ref_id, ref_name, region):
        igv_notebook.init()

        if ref_fname_or_genome_name.endswith(
            (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")
        ):
            self.add_genome_from_fasta(
                ref_fname_or_genome_name, ref_id, ref_name, region
            )
        else:
            self.add_genome(ref_fname_or_genome_name, region)

    def add_genome_from_fasta(
        self,
        ref_fname: str,
        ref_id: Optional[str] = None,
        ref_name: Optional[str] = None,
        region: Optional[str] = None,
    ):
        if ref_id is not None and ref_name is None:
            ref_name = ref_id
        elif ref_id is None and ref_name is not None:
            ref_id = ref_name
        else:
            logger.error("One of `ref_id` and `ref_name` must be not None.")
        self.b = igv_notebook.Browser(
            {
                "reference": {
                    "id": ref_id,
                    "name": ref_name,
                    "fastaURL": ref_fname,
                    "indexURL": f"{ref_fname}.fai",
                },
                "locus": region,
            }
        )
        return self

    def add_genome(self, genome_name: str, region: Optional[str] = None):
        """`genome_name` is e.g. 'hg38'"""
        self.b = igv_notebook.Browser({"genome": genome_name, "locus": region})
        return self

    def add_bam(
        self,
        bam_fname: str,
        name: Optional[str] = "bam",
        height: Optional[int] = 300,
        auto_height: bool = True,
    ):
        self.b.load_track(
            {
                "name": name,
                "path": bam_fname,
                "indexPath": f"{bam_fname}.bai",
                "format": "bam",
                "type": "alignment",
                "height": height,
                "autoHeight": auto_height,
            }
        )
        return self

    def add_vcfgz(
        self,
        vcfgz_fname: str,
        name: Optional[str] = "vcf",
        height: Optional[int] = 50,
        auto_height: bool = True,
    ):
        self.b.load_track(
            {
                "name": name,
                "path": vcfgz_fname,
                "indexPath": f"{vcfgz_fname}.tbi",
                "format": "vcf",
                "type": "variant",
                "height": height,
                "autoHeight": auto_height,
            }
        )
        return self

    def add_gff3gz(
        self,
        gff3gz_fname: str,
        name: Optional[str] = "gff",
        height: Optional[int] = 50,
        auto_height: bool = True,
    ):
        self.b.load_track(
            {
                "name": name,
                "path": gff3gz_fname,
                "indexPath": f"{gff3gz_fname}.tbi",
                "format": "gff3",
                "type": "annotation",
                "height": height,
                "autoHeight": auto_height,
            }
        )
        return self

    def add_bed(
        self,
        bed_fname: str,
        name: Optional[str] = None,
        height: Optional[int] = 50,
        auto_height: bool = True,
    ):
        """bed, bb, bed.gz, gtf are accepted."""
        ext = bed_fname.split(".")[-1]
        is_gz = False
        if ext == "gz":
            ext = bed_fname.split(".")[-2]
            is_gz = True
        if name is None:
            name = ext
        if not is_gz:
            self.b.load_track(
                {
                    "name": name,
                    "path": bed_fname,
                    "format": ext,
                    "type": "annotation",
                    "height": height,
                    "autoHeight": auto_height,
                }
            )
        else:
            self.b.load_track(
                {
                    "name": name,
                    "path": bed_fname,
                    "indexPath": f"{bed_fname}.tbi",
                    "format": ext,
                    "type": "annotation",
                    "height": height,
                    "autoHeight": auto_height,
                }
            )
        return self

    def add_wig(
        self,
        wig_fname: str,
        name: Optional[str] = "annot",
        height: Optional[int] = 50,
        auto_height: bool = True,
    ):
        """wig, bigwig, bedgraph are accepted."""
        self.b.load_track(
            {
                "name": name,
                "path": wig_fname,
                "type": "wig",
                "height": height,
                "autoHeight": auto_height,
            }
        )
        return self
