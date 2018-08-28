from logzero import logger

edlib_mode = {"global": "NW", "glocal": "HW", "local": "SHW"}


def cigar_to_len(cigar):
    alignment_len = 0
    d = ''
    for c in cigar:
        if c == '=' or c == 'D' or c == 'I' or c == 'X':
            alignment_len += int(d)
            d = ''
        else:
            d += c
    return alignment_len


def revcomp_cigar(cigar):
    ret = ""
    while len(cigar) > 0:
        cut_point = 0
        while cigar[cut_point] not in set(['=', 'X', 'I', 'D']):
            cut_point += 1
        ret = cigar[0:cut_point + 1] + ret
        cigar = cigar[cut_point + 1:]
    return ret


def run_edlib(query,
              target,
              mode,
              only_diff=False,
              revcomp=False,
              revcomp_pruning_diff_th=0.3,   # at most this diff if true strand
              strand_prior=None,   # 0 or 1
              cyclic=False,
              return_seq=False,   # in glocal alignment
              return_seq_diff_th=1.0,   # remove distant sequences
              return_cigar=False,
              swap_cigar=False):
    """
    Perform 1 vs 1 pairwise sequence alignment using edlib.
    <mode> = {"global", "glocal", "local"}
    In "glocal" mode, boundary-fixed <query> is mapped to <target>.

    Originally, "<query> == cigar(<target>[<start>:<end>])" holds.
    If <swap_cigar> is True, "cigar(<query>) == <target>[<start>:<end>]" does.
    This is effective when you actually want to consider <query> as target
    (i.e. reference) in the SAM file format.

    NOTE: local alignment using edlib is not recommended.
    """

    # TODO: implement <revcomp_query> (or <multi_target>?) and  <count_bad_align>

    assert mode in ["global", "glocal", "local"], f"Invalid mode name: {mode}"
    assert query != "" and target != "", "Input sequence must not be empty!"

    import edlib
    if revcomp:
        from .seq import revcomp

    # NOTE: glocal cyclic can return a sequence shorter than the original one
    global_cyclic = False
    if cyclic:
        target = target * 2
        if mode == "global":
            global_cyclic = True
            mode = "glocal"

    mode = edlib_mode[mode]

    strand = 0
    if not revcomp:
        align = edlib.align(query, target, mode=mode, task="path")
        diff = align["editDistance"] / cigar_to_len(align["cigar"])
    else:
        if strand_prior is None or strand_prior == 0:
            align_f = edlib.align(query, target, mode=mode, task="path")
            diff_f = align_f["editDistance"] / cigar_to_len(align_f["cigar"])
            if diff_f < revcomp_pruning_diff_th:
                align = align_f
                diff = diff_f
            else:
                target_rc = revcomp(target)
                align_rc = edlib.align(query, target_rc, mode=mode, task="path")
                diff_rc = align_rc["editDistance"] / cigar_to_len(align_rc["cigar"])
                if diff_f <= diff_rc:
                    align = align_f
                    diff = diff_f
                else:
                    strand = 1
                    align = align_rc
                    diff = diff_rc
                    target = target_rc
        else:
            target_rc = revcomp(target)
            align_rc = edlib.align(query, target_rc, mode=mode, task="path")
            diff_rc = align_rc["editDistance"] / cigar_to_len(align_rc["cigar"])
            if diff_rc < revcomp_pruning_diff_th:
                strand = 1
                align = align_rc
                diff = diff_rc
                target = target_rc
            else:
                align_f = edlib.align(query, target, mode=mode, task="path")
                diff_f = align_f["editDistance"] / cigar_to_len(align_f["cigar"])
                if diff_f <= diff_rc:
                    align = align_f
                    diff = diff_f
                else:
                    strand = 1
                    align = align_rc
                    diff = diff_rc
                    target = target_rc

    # NOTE: coordinates are on <target>
    start, end = align["locations"][0]
    end += 1   # for compatibility with python's slice

    flag_cycle = False
    if cyclic:
        half_pos = len(target) // 2
        if end >= half_pos:
            end -= half_pos
            flag_cycle = not flag_cycle
            if start >= half_pos:
                start -= half_pos
                flag_cycle = not flag_cycle
    
        if global_cyclic:
            if flag_cycle:
                n_dist = start - end   # >0: ins, <0: del at boundary
                end = start
            else:   # adjust so that the alignment is same as global mode
                n_dist = start + half_pos - end
                start = 0
                end = half_pos
            diff = (align["editDistance"] + n_dist) / (cigar_to_len(align["cigar"]) + n_dist)
            # TODO: adjust cigar at bounary (deletion case is annoying)
            
        target = target[:half_pos]

    if only_diff:
        return diff

    if return_seq:
        if diff > return_seq_diff_th:
            seq = None
        else:
            if flag_cycle:
                seq = target[start:] + target[:end]
            else:
                seq = target[start:end]

    if strand == 1:
        start_tmp = start
        start = len(target) - end
        end = len(target) - start_tmp

    ret = {"start": start, "end": end, "diff": diff, "strand": strand}

    if return_seq:
        ret["seq"] = seq

    if return_cigar:
        cigar = (align["cigar"] if not swap_cigar
                 else align["cigar"].replace("I", "?").replace("D", "I").replace("?", "D"))

        # NOTE: "cigar(revcomp(target[start:end])) == query" when strand == 1
        ret["cigar"] = cigar

        # NOTE: if you prefer "revcomp(cigar(target[start:end])) == query", then use below
        #ret["cigar"] = cigar if strand == 0 else revcomp_cigar(cigar)

    return ret


"""
def run_consed(in_seqs,
               only_consensus=True,
               out_prefix="out",
               variant_vector=False,
               variant_graph=False,
               variant_fraction=0.3,
               display_width=100000,
               tmp_dir="tmp",
               n_iteration=1,   # times to iteratively run Consed
               parallel=False):   # avoid file name collision
    # TODO: implement second stage of consed

    in_fname = f"{tmp_dir}/consed.seqs"
    if parallel:
        in_fname += f".{os.getpid()}"   # TODO: implement parallel mode for not only_consensu mode
    with open(in_fname, 'w') as f:
        f.write('\n'.join(in_seqs) + '\n')

    if only_consensus:
        command = f"consed {in_fname}"
    else:
        vv = '-V' if variant_vector else ''
        vg = f"-G{out_prefix}" if variant_graph else ''
        vf = f"-t{variant_fraction}" if variant_vector or variant_graph else ''
        command = f"consed {vv} {vg} {vf} -w{display_width} {in_fname}"

    consed_out = run_command(command)   # TODO: maybe here should also be separated into each mode

    if only_consensus:
        consed_out = consed_out.strip().split('\n')
        start_index = 0
        while consed_out[start_index][0] not in set(['a', 'c', 'g', 't']):
            start_index += 1
            if start_index == len(consed_out):
                logger.error("No consensus sequence in the output")
                return ''
        if start_index > 0:
            logger.debug(f"{consed_out[:start_index]}")
        return ''.join(consed_out[start_index:])
    else:
        with open(f"{out_prefix}.consed", 'w') as f:
            f.write(f"{consed_out}")

        # Extract consensus sequence
        command = f"awk 'BEGIN {{seq = \"\"}} $0 == \"\" {{exit}} {{seq = seq $0}} END {{print seq}}' {out_prefix}.consed"
        consensus_seq = run_command(command).strip()
        with open(f"{out_prefix}.consensus", 'w') as f:
            f.write(f">consensus/0/0_{len(consensus_seq)}\n{consensus_seq}\n")

        # Extract variant matrix
        command = f"sed -e '0,/Variant/d' {out_prefix}.consed | sed -e '1,2d' | awk 'BEGIN {{i = 1}} {{if ($0 == \"\") {{i = 1}} else {{seq[i] = seq[i] $NF; i++}}}} END {{for (i = 1; i <= length(seq); i++) print seq[i]}}'"
        variant_matrix = run_command(command)
        with open(f"{out_prefix}.V", 'w') as f:
            f.write(f"{variant_matrix}")
"""
