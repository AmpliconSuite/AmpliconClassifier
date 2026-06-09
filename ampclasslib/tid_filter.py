"""
Detect tandem inverted duplications (TIDs) misclassified as ecDNA.

Criteria:
  C1: Exactly one -- and one ++ discordant edge on the same chromosome.
  C2: No other discordant edges above noise thresholds where both endpoints are amplified (CN >= C2_AMP_CN_FLOOR).
  C3: Max inner-segment CN <= CN_RATIO_MAX * outer-flank CN.
  C4: All circular cycles have segments contained within the TID region.
  C5: The -- and ++ spans do not significantly overlap.

Primary entry point: check_tid(seq_edges, sv_edges, cycleList, segSeqD)
  seq_edges, sv_edges — as returned by parse_graph_edges_raw()
  cycleList, segSeqD  — as returned by parseCycle()
"""

MIN_SV_SIZE             = 10_000
CN_RATIO_MAX            = 2.5
CN_RATIO_MAX_TIGHT_FB   = 4.0     # relaxed C3 cap when one foldback is tight (< TIGHT_FB_SIZE)
TIGHT_FB_SIZE           = 50_000  # foldback span below this qualifies as tight
MIN_DISC_CN             = 0.5
MIN_INTERCHROM_CN       = 0.2
MAX_EDGE_OVERLAP        = 0.10
MIN_SMALL_EDGE          = 5_000
C2_AMP_CN_FLOOR         = 4.0     # SV endpoint with CN below this is escaping the amplicon; ignore for C2
C2_MIN_BP_CN_DELTA      = 2.0     # if neither endpoint of an other-SV shows a CN step >= this, treat as noise
C4_MIN_NONINV_CYCLE_SIZE = 50_000  # non-inversion circular cycles smaller than this are ignored in C4
MAX_FOLDBACK_SPAN       = 1_000_000  # reject if either selected foldback spans more than this


def check_tid(seq_edges, sv_edges, cycleList, segSeqD,
              cn_ratio_max=CN_RATIO_MAX,
              cn_ratio_max_tight_fb=CN_RATIO_MAX_TIGHT_FB,
              tight_fb_size=TIGHT_FB_SIZE,
              max_foldback_span=MAX_FOLDBACK_SPAN):
    """
    Returns a dict with 'pass': True on success, or 'fail': <reason> on failure,
    plus diagnostic fields (ratio, max_inner, background, tid_span_kb, etc.).
    """
    # Translate sv_edges dicts to flat tuples (c1, p1, s1, c2, p2, s2, cn, _)
    disc = [
        (e['chrom1'], e['pos1'], e['strand1'],
         e['chrom2'], e['pos2'], e['strand2'],
         e['cn'], 0)
        for e in sv_edges
    ]

    mm_raw, pp_raw, other = [], [], []
    for d in disc:
        c1, p1, s1, c2, p2, s2, cn, _ = d
        if c1 != c2:
            if cn >= MIN_INTERCHROM_CN:
                other.append(d)
            continue
        if cn < MIN_DISC_CN:
            continue
        span = abs(p2 - p1)
        if s1 == '-' and s2 == '-':
            mm_raw.append(d)
        elif s1 == '+' and s2 == '+':
            pp_raw.append(d)
        elif span > MIN_SV_SIZE:
            other.append(d)

    def drop_small_secondary(edges):
        if len(edges) <= 1:
            return edges
        large = [e for e in edges if abs(e[1] - e[4]) >= MIN_SMALL_EDGE]
        return large if large else edges

    # Build the set of breakpoint vertex pairs actually traversed in cycles.
    def _cycle_breakpoints():
        used = set()
        for cycle in cycleList:
            n = len(cycle)
            for i in range(n):
                si, sj = cycle[i], cycle[(i + 1) % n]
                si_id, sj_id = abs(si), abs(sj)
                if si_id == 0 or sj_id == 0:
                    continue
                if si_id not in segSeqD or sj_id not in segSeqD:
                    continue
                si_c, si_s, si_e = segSeqD[si_id]
                sj_c, sj_s, sj_e = segSeqD[sj_id]
                si_dir = '+' if si > 0 else '-'
                sj_dir = '+' if sj > 0 else '-'
                exit_p = si_e if si_dir == '+' else si_s
                exit_s = '+' if si_dir == '+' else '-'
                entry_p = sj_s if sj_dir == '+' else sj_e
                entry_s = '-' if sj_dir == '+' else '+'
                used.add((si_c, exit_p, exit_s, sj_c, entry_p, entry_s))
        return used

    cycle_bps = _cycle_breakpoints()

    def _in_cycles(d):
        k1 = (d[0], d[1], d[2], d[3], d[4], d[5])
        k2 = (d[3], d[4], d[5], d[0], d[1], d[2])
        return k1 in cycle_bps or k2 in cycle_bps

    def _cycle_filter(edges):
        if not cycle_bps:
            return edges
        filtered = [e for e in edges if _in_cycles(e)]
        return filtered if filtered else edges

    def _bp_cn_quick(chrom, pos, strand):
        for e in seq_edges:
            if e['chrom'] != chrom:
                continue
            if strand == '+' and e['end'] == pos:
                return e['cn']
            if strand == '-' and e['start'] == pos:
                return e['cn']
        return None

    def _amp_filter(edges):
        """Drop foldback edges where both endpoints are outside the amplicon (CN < C2_AMP_CN_FLOOR)."""
        def in_amp(d):
            cn1 = _bp_cn_quick(d[0], d[1], d[2])
            cn2 = _bp_cn_quick(d[3], d[4], d[5])
            return (cn1 is not None and cn1 >= C2_AMP_CN_FLOOR) or \
                   (cn2 is not None and cn2 >= C2_AMP_CN_FLOOR)
        filtered = [e for e in edges if in_amp(e)]
        return filtered if filtered else edges

    mm = drop_small_secondary(_amp_filter(_cycle_filter(mm_raw)))
    pp = drop_small_secondary(_amp_filter(_cycle_filter(pp_raw)))

    def fmt(d):
        c1, p1, s1, c2, p2, s2, cn, _ = d
        span = abs(p2 - p1) if c1 == c2 else float('inf')
        return "{c1}:{p1}{s1}->{c2}:{p2}{s2} CN={cn:.2f} span={span:.1f}kb".format(
            c1=c1, p1=p1, s1=s1, c2=c2, p2=p2, s2=s2, cn=cn, span=span / 1000)

    # Helper functions needed for both pre-C1 foldback filtering and C2.
    def _bp_cn(chrom, pos, strand):
        for e in seq_edges:
            if e['chrom'] != chrom:
                continue
            if strand == '+' and e['end'] == pos:
                return e['cn']
            if strand == '-' and e['start'] == pos:
                return e['cn']
        return None

    def _opposite_bp_cn(chrom, pos, strand):
        for e in seq_edges:
            if e['chrom'] != chrom:
                continue
            if strand == '+' and e['start'] == pos + 1:
                return e['cn']
            if strand == '-' and e['end'] == pos - 1:
                return e['cn']
        return None

    def _bp_cn_delta(chrom, pos, strand):
        cn1 = _bp_cn(chrom, pos, strand)
        cn2 = _opposite_bp_cn(chrom, pos, strand)
        if cn1 is None or cn2 is None:
            return float('inf')
        return abs(cn1 - cn2)

    def _in_amplicon(chrom, pos, strand):
        cn = _bp_cn(chrom, pos, strand)
        return cn is not None and cn >= C2_AMP_CN_FLOOR

    def _creates_cn_change(d):
        return (_bp_cn_delta(d[0], d[1], d[2]) >= C2_MIN_BP_CN_DELTA or
                _bp_cn_delta(d[3], d[4], d[5]) >= C2_MIN_BP_CN_DELTA)

    # When multiple -- or ++ edges exist, drop CN-neutral ones (secondary structural elements).
    # If all edges are CN-neutral (filter would empty the set), fall back to the highest-CN edge.
    if len(mm) > 1:
        filtered = [d for d in mm if _creates_cn_change(d)]
        mm = filtered if filtered else [max(mm, key=lambda d: d[6])]
    if len(pp) > 1:
        filtered = [d for d in pp if _creates_cn_change(d)]
        pp = filtered if filtered else [max(pp, key=lambda d: d[6])]

    info = {
        'n_mm': len(mm), 'n_pp': len(pp), 'n_other': len(other),
        'mm_str': [fmt(d) for d in mm],
        'pp_str': [fmt(d) for d in pp],
        'other_str': [fmt(d) for d in other],
    }

    # C1
    if len(mm) != 1 or len(pp) != 1:
        info['fail'] = "C1: {} -- edge(s), {} ++ edge(s) (need 1 each)".format(len(mm), len(pp))
        return info

    # C2: SVs where either endpoint lands in unamplified sequence (CN < C2_AMP_CN_FLOOR)
    # are escaping the amplicon and don't affect the internal CN structure — ignore them.
    other_internal = [d for d in other
                      if _in_amplicon(d[0], d[1], d[2]) and _in_amplicon(d[3], d[4], d[5])
                      and _creates_cn_change(d)]
    info['other_internal_str'] = [fmt(d) for d in other_internal]

    if other_internal:
        info['fail'] = "C2: {} internal large SV(s): {}".format(
            len(other_internal), [fmt(d) for d in other_internal])
        return info

    mm_e, pp_e = mm[0], pp[0]
    chrom = mm_e[0]
    if pp_e[0] != chrom:
        info['fail'] = "C1: -- and ++ on different chromosomes"
        return info

    # C1b: neither foldback may span more than max_foldback_span.
    # Multi-Mb foldbacks are not biologically consistent with a simple TID.
    mm_span = abs(mm_e[4] - mm_e[1])
    pp_span = abs(pp_e[4] - pp_e[1])
    if mm_span > max_foldback_span or pp_span > max_foldback_span:
        info['fail'] = "C1b: foldback span exceeds {}kb (-- {:.0f}kb, ++ {:.0f}kb)".format(
            max_foldback_span // 1000, mm_span / 1000, pp_span / 1000)
        return info

    # C5: -- and ++ spans must not significantly overlap
    def edge_overlap_frac(e1, e2):
        a, b = sorted([e1[1], e1[4]])
        c, d = sorted([e2[1], e2[4]])
        overlap = max(0, min(b, d) - max(a, c))
        shorter = min(b - a, d - c)
        return overlap / shorter if shorter > 0 else 0.0

    min_frac = min(edge_overlap_frac(me, pe) for me in mm_raw for pe in pp_raw)
    if min_frac > MAX_EDGE_OVERLAP:
        info['fail'] = "C5: edges cross (min overlap fraction {:.2f} > {})".format(
            min_frac, MAX_EDGE_OVERLAP)
        return info

    # Outermost breakpoints across all raw same-orientation edges
    all_bps = sorted(
        [(d[1], d[2]) for d in mm_raw + pp_raw] +
        [(d[4], d[5]) for d in mm_raw + pp_raw],
        key=lambda x: x[0]
    )
    lpos, lstrand = all_bps[0]
    rpos, rstrand = all_bps[-1]

    # C3: outer-flank CNs.
    # Use the minimum CN across ALL segments on each side of the TID boundary — not just the
    # immediately adjacent one. A small sub-threshold SV between the foldback endpoint and the
    # true diploid flank can leave a still-amplified adjacent segment; taking the min finds the
    # genuine diploid flank further out.
    left_boundary  = lpos - 1 if lstrand == '-' else lpos
    right_boundary = rpos + 1 if rstrand == '+' else rpos

    left_segs  = [e['cn'] for e in seq_edges if e['chrom'] == chrom and e['end']   <= left_boundary]
    right_segs = [e['cn'] for e in seq_edges if e['chrom'] == chrom and e['start'] >= right_boundary]

    left_cn  = min(left_segs)  if left_segs  else None
    right_cn = min(right_segs) if right_segs else None

    info.update({'left_cn': left_cn, 'right_cn': right_cn,
                 'tid_span_kb': (rpos - lpos) / 1000})

    if left_cn is None or right_cn is None:
        info['fail'] = "C3: flank lookup failed (left={}, right={})".format(left_cn, right_cn)
        return info

    inner = [e for e in seq_edges
             if e['chrom'] == chrom and e['start'] >= lpos and e['end'] <= rpos]
    if not inner:
        info['fail'] = "C3: no inner segments"
        return info

    total_size = sum(e['end'] - e['start'] for e in inner)
    inner_cn   = sum(e['cn'] * (e['end'] - e['start']) for e in inner) / total_size
    background = max(left_cn, right_cn)
    ratio      = inner_cn / background

    tight_fb = min(mm_span, pp_span) < tight_fb_size
    ratio_max = cn_ratio_max_tight_fb if tight_fb else cn_ratio_max

    info.update({'inner_cn': inner_cn, 'background': background, 'ratio': ratio,
                 'tight_fb': tight_fb, 'ratio_max': ratio_max})

    if ratio > ratio_max:
        info['fail'] = "C3: ratio {:.2f}x > {:.1f}x (inner={:.2f}, bg={:.2f}{})".format(
            ratio, ratio_max, inner_cn, background,
            ', tight_fb' if tight_fb else '')
        return info

    # C3b: the sequence segment immediately inside each foldback endpoint must be amplified.
    # For a real TID the amplified region starts at the foldback boundaries; if any inner face
    # has CN < 0.5 * inner_cn the foldback is not sitting at the edge of the amplified region.
    def _inner_face_cn(pos, strand):
        """CN of the seq segment on the inner (TID-region) side of a foldback endpoint."""
        for e in seq_edges:
            if e['chrom'] != chrom:
                continue
            if strand == '-' and e['start'] == pos:   # edge exits left → inner = right
                return e['cn']
            if strand == '+' and e['end'] == pos:     # edge exits right → inner = left
                return e['cn']
        return None

    face_threshold = inner_cn * 0.5
    face_cns = {
        'mm_p1': _inner_face_cn(mm_e[1], mm_e[2]),
        'mm_p2': _inner_face_cn(mm_e[4], mm_e[5]),
        'pp_p1': _inner_face_cn(pp_e[1], pp_e[2]),
        'pp_p2': _inner_face_cn(pp_e[4], pp_e[5]),
    }
    info['face_cns'] = face_cns
    low_faces = {k: v for k, v in face_cns.items() if v is not None and v < face_threshold}
    if low_faces:
        info['fail'] = (
            "C3b: foldback inner face(s) below {:.2f} (0.5x inner_cn={:.2f}): {}".format(
                face_threshold, inner_cn,
                ", ".join("{}={:.2f}".format(k, v) for k, v in low_faces.items()))
        )
        return info

    # C4: all circular cycles must be contained within the TID region.
    # Non-inversion cycles (all segments same orientation) under C4_MIN_NONINV_CYCLE_SIZE
    # are byproduct circlets and are ignored.
    ext_pos = ([d[1] for d in mm_raw + pp_raw] + [d[4] for d in mm_raw + pp_raw])
    ext_lpos, ext_rpos = min(ext_pos), max(ext_pos)

    def _cycle_size(cycle):
        return sum(segSeqD[abs(s)][2] - segSeqD[abs(s)][1]
                   for s in cycle if s != 0 and abs(s) in segSeqD)

    def _is_inversion_cycle(cycle):
        dirs = set(1 if s > 0 else -1 for s in cycle if s != 0)
        return len(dirs) > 1

    n_circ = sum(1 for cycle in cycleList if 0 not in cycle)
    bad = []
    for cycle in cycleList:
        if 0 in cycle:
            continue
        if not _is_inversion_cycle(cycle) and _cycle_size(cycle) < C4_MIN_NONINV_CYCLE_SIZE:
            continue
        for seg_id in cycle:
            sid = abs(seg_id)
            if sid == 0 or sid not in segSeqD:
                continue
            sc, ss, se = segSeqD[sid]
            if sc != chrom or ss < ext_lpos or se > ext_rpos:
                bad.append(sid)
                break

    info['n_circular'] = n_circ
    if bad:
        info['fail'] = "C4: {} circular cycle(s) outside TID region".format(len(bad))
        return info

    info['pass'] = True
    return info
