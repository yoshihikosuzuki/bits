paf_fname = "/home/yoshihiko_s/work/medaka/Hd-rR/20-scaffolds/v3.0.2/40-patch-misasm/chr1_37643270_37703026/scaffolds.ont.wm.chr1_37643270_37703026.backbone.others.round5.paf"

data = list(
    map(
        lambda line: line.split("\t"),
        bu.run_command(f"cat {paf_fname}").strip().split("\n"),
    )
)


qn, ql, qb, qe, s, tn, tl, tb, te = data[0][:9]

pl.show(
    pl.lines(
        [
            (int(qb), int(tb), int(qe), int(te))
            for qn, ql, qb, qe, s, tn, tl, tb, te in map(lambda line: line[:9], data)
        ],
        use_webgl=False,
    ),
    pl.layout(
        width=600 + 50,
        height=600 / int(ql) * int(tl) + 90,
        x_range=(0, ql),
        y_range=(0, tl),
        x_zeroline=True,
        y_zeroline=True,
        x_title="Query",
        y_title="Target",
        x_grid=False,
        y_grid=False,
        x_bounding_line=True,
        x_mirror=True,
        y_bounding_line=True,
        y_mirror=True,
        anchor_axes=True,
    ),
)
