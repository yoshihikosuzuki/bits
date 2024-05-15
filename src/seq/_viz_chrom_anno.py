# NOTE: Visualization for per-chromosome annotation data

import plotly_light as pl


def detail_depth_to_trace(data, mean_cov, max_cov, chrom_len, title):
    x = [r.b for r in data]
    return pl.figure(
        [
            pl.lines(
                (0, mean_cov, chrom_len, mean_cov),
                width=2,
                col="gray",
                name="global mean",
                use_webgl=False,
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.min for r in data],
                mode="lines",
                col=pl.colors["red"],
                name="min",
                use_webgl=False,
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.max for r in data],
                mode="lines",
                col=pl.colors["yellow"],
                name="max",
                use_webgl=False,
                show_legend=True,
            ),
            pl.scatter(
                x,
                [r.med for r in data],
                mode="lines",
                col=pl.colors["blue"],
                name="median",
                line_width=2,
                use_webgl=False,
                show_legend=True,
            ),
        ],
        pl.layout(
            title=title,
            y_range=(0, max_cov),
            x_bounding_line=True,
            y_bounding_line=True,
            x_mirror=True,
            y_mirror=True,
            x_grid=True,
            y_grid=True,
        ),
    )


def mean_depth_to_trace(data, mean_cov, max_cov, chrom_len, title):
    return pl.figure(
        [
            pl.lines(
                (0, mean_cov, chrom_len, mean_cov), width=1, col="gray", use_webgl=False
            ),
            pl.scatter(
                [r.b for r in data],
                [r.cov for r in data],
                mode="lines",
                col="gray",
                line_width=1,
                use_webgl=False,
            ),
        ],
        pl.layout(
            title=title,
            y_range=(0, max_cov),
            x_bounding_line=True,
            y_bounding_line=True,
            x_mirror=True,
            y_mirror=True,
            x_grid=True,
            y_grid=True,
        ),
    )


def bed_to_trace(data, pad, chrom, chrom_len, name, col):
    data_filt = list(
        filter(lambda x: x.chr == chrom and x.e > 0 and x.b < chrom_len, data)
    )
    if len(data_filt) > 0:
        coords, texts = map(
            lambda x: list(x),
            zip(
                *[
                    ((max(x.b - pad, 0), name, x.e + pad, name), f"{x.chr}:{x.b}-{x.e}")
                    for x in data_filt
                ]
            ),
        )
        return pl.lines(
            coords, text=texts, width=8, col=col, name=name, show_legend=False
        )
    else:
        return pl.lines(
            [(0, name, 0, name)], width=0, col=col, name=name, show_legend=False
        )


@dataclass
class Annot:
    name: str
    data: List[bs.BedRecord]
    col: str
    pad: float = 0.0005


@dataclass
class AsmPlotter:
    def __post_init__(self):
        root = f"{strain}/20-scaffolds/{version}"
        self.cov = {
            "ont": bu.load_pickle(f"jupyter_pickle/{strain}.{version}.ont.cov.bin.pkl"),
            "hifi male": bs.load_bed(
                f"{root}/04-winnowmap-hifi/scaffolds.hifi.male.m10k.wm.sorted.bam.regions.bedgraph",
                attrs=[("cov", float)],
            ),
            "hifi female": bs.load_bed(
                f"{root}/04-winnowmap-hifi/scaffolds.hifi.female.m10k.wm.sorted.bam.regions.bedgraph",
                attrs=[("cov", float)],
            ),
        }
        self.annot = [
            Annot(
                "sequence gap",
                bs.load_bed(f"{root}/07-asset/scaffolds.gaps.bed"),
                "black",
                0.001,
            ),
            Annot(
                "centromere",
                bs.load_bed(f"{root}/30-trf/scaffolds.cen.bed"),
                pl.colors["darkyellow"],
            ),
            Annot(
                "telomere",
                bs.load_bed(f"{root}/30-trf/scaffolds.tel.bed"),
                pl.colors["darkyellow"],
            ),
            Annot(
                "TR",
                bs.load_bed(f"{root}/30-trf/scaffolds.tr.bed"),
                pl.colors["darkyellow"],
                0,
            ),
            Annot(
                "45S rDNA",
                bs.load_bed(
                    f"{root}/34-rDNA-45S/45S.scaffolds.minimap2.asm20.filtered.bed"
                ),
                pl.colors["darkyellow"],
            ),
            Annot(
                "5S rDNA",
                bs.load_bed(
                    f"{root}/34-rDNA-5S/5S.scaffolds.minimap2.asm20.filtered.bed"
                ),
                pl.colors["darkyellow"],
            ),
            Annot(
                "teratorn",
                bs.load_bed(
                    f"{root}/31-teratorn/teratorn.scaffolds.minimap2.asm20.bed"
                ),
                pl.colors["darkyellow"],
            ),
            Annot(
                "tol2",
                bs.load_bed(f"{root}/32-tol2/tol2.scaffolds.minimap2.asm20.bed"),
                pl.colors["blue"],
            ),
            Annot(
                "HiFi SNV",
                bs.load_bed(
                    f"{root}/05-deepvariant-hifi/scaffolds.hifi.wm.dv.pass.snp.bed"
                ),
                pl.colors["red"],
            ),
            #                      Annot("ONT SNV", bs.load_bed(f"{root}/05-deepvariant-ont/scaffolds.ont.wm.dv.pass.snp.bed"), pl.colors["red"]),
            Annot(
                "ONT variants (GQ>25)",
                bs.load_bed(
                    f"{root}/05-deepvariant-ont/scaffolds.ont.wm.dv.pass.GQ25.bed"
                ),
                pl.colors["red"],
            ),
            #                      Annot("ONT patch region", bs.load_bed(f"{root}/40-patch-misasm/patch_regions.bed"), pl.colors["purple"], 0.0001),
            #                      Annot("ONT patch region (variant)", bs.load_bed(f"{root}/40-patch-misasm/regions.variant.bed"), pl.colors["lightpurple"], 0.0001),
            #                      Annot("ONT patch region (gap)", bs.load_bed(f"{root}/40-patch-misasm/regions.gap.bed"), pl.colors["lightpurple"], 0.0001),
            #                      Annot("ONT patch region (cov)", bs.load_bed(f"{root}/40-patch-misasm/regions.ab_cov.ont.not_gap.not_variant.bed"), pl.colors["lightpurple"], 0.0001),
        ]

    def plot(
        self,
        chrom,
        ont_max_cov=100,
        hifi_max_cov=80,
        annot_width_ratio=1.0,
        do_output=False,
    ):
        chrom_len = scafs_by_name[chrom].length
        fig_d_ont = detail_depth_to_trace(
            self.cov["ont"][chrom],
            mean_cov_ont[strain],
            ont_max_cov,
            chrom_len,
            "ONT UL coverage",
        )
        fig_d_hifi_male = mean_depth_to_trace(
            list(filter(lambda x: x.chr == chrom, self.cov["hifi male"])),
            mean_cov_hifi[strain] / 2,
            hifi_max_cov,
            chrom_len,
            "HiFi mean coverage (male)",
        )
        fig_d_hifi_female = mean_depth_to_trace(
            list(filter(lambda x: x.chr == chrom, self.cov["hifi female"])),
            mean_cov_hifi[strain] / 2,
            hifi_max_cov,
            chrom_len,
            "HiFi mean coverage (female)",
        )
        fig_a = pl.figure(
            [
                bed_to_trace(
                    a.data,
                    chrom_len * a.pad * annot_width_ratio,
                    chrom,
                    chrom_len,
                    a.name,
                    a.col,
                )
                for a in self.annot
            ],
            pl.layout(y_reversed=True, y_category=True, y_tickfontsize=12, y_grid=True),
        )
        if do_output:
            dname = f"jupyter_plot/coverage/{strain}/{version}"
            bu.run_command(f"mkdir -p {dname}")
        pl.show_mult(
            [fig_d_ont, fig_d_hifi_female, fig_d_hifi_male, fig_a],
            n_col=1,
            row_heights=[1, 0.5, 0.5, 1.25],
            vertical_spacing=0.05,
            shared_xaxes=True,
            layout=pl.layout(
                width=1300, height=800, title=f"{chrom}", x_range=(0, chrom_len)
            ),
            out_html=f"{dname}/{strain}.{chrom}.html" if do_output else None,
            out_image=f"{dname}/{strain}.{chrom}.pdf" if do_output else None,
            return_fig=do_output,
        )
