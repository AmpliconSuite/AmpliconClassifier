"""
Detect tandem inverted duplications (TIDs) misclassified as ecDNA.

Criteria:
  C1: Exactly one -- and one ++ discordant edge on the same chromosome.
  C2: No other discordant edges above noise thresholds.
  C3: Max inner-segment CN <= CN_RATIO_MAX * outer-flank CN.
  C4: All circular cycles have segments contained within the TID region.
  C5: The -- and ++ spans do not significantly overlap.

Primary entry point: check_tid(seq_edges, sv_edges, cycleList, segSeqD)
  seq_edges, sv_edges — as returned by parse_graph_edges_raw()
  cycleList, segSeqD  — as returned by parseCycle()
"""

MIN_SV_SIZE       = 10_000
CN_RATIO_MAX      = 2.5
MIN_DISC_CN       = 0.5
MIN_INTERCHROM_CN = 0.2
MAX_EDGE_OVERLAP  = 0.10
MIN_SMALL_EDGE    = 5_000


def check_tid(seq_edges, sv_edges, cycleList, segSeqD):
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

    mm = drop_small_secondary(mm_raw)
    pp = drop_small_secondary(pp_raw)

    def fmt(d):
        c1, p1, s1, c2, p2, s2, cn, _ = d
        span = abs(p2 - p1) if c1 == c2 else float('inf')
        return "{c1}:{p1}{s1}->{c2}:{p2}{s2} CN={cn:.2f} span={span:.1f}kb".format(
            c1=c1, p1=p1, s1=s1, c2=c2, p2=p2, s2=s2, cn=cn, span=span / 1000)

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
    # C2
    if other:
        info['fail'] = "C2: {} other large SV(s): {}".format(len(other), [fmt(d) for d in other])
        return info

    mm_e, pp_e = mm[0], pp[0]
    chrom = mm_e[0]
    if pp_e[0] != chrom:
        info['fail'] = "C1: -- and ++ on different chromosomes"
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

    # C3: outer-flank CNs
    def seg_by_end(chrom, end):
        for e in seq_edges:
            if e['chrom'] == chrom and e['end'] == end:
                return e['cn']
        return None

    def seg_by_start(chrom, start):
        for e in seq_edges:
            if e['chrom'] == chrom and e['start'] == start:
                return e['cn']
        return None

    left_cn  = seg_by_end(chrom, lpos - 1) if lstrand == '-' else seg_by_end(chrom, lpos)
    right_cn = seg_by_start(chrom, rpos + 1) if rstrand == '+' else seg_by_start(chrom, rpos)

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

    max_inner  = max(e['cn'] for e in inner)
    background = max(left_cn, right_cn)
    ratio      = max_inner / background

    info.update({'max_inner': max_inner, 'background': background, 'ratio': ratio})

    if ratio > CN_RATIO_MAX:
        info['fail'] = "C3: ratio {:.2f}x > {:.1f}x (inner={:.2f}, bg={:.2f})".format(
            ratio, CN_RATIO_MAX, max_inner, background)
        return info

    # C4: all circular cycles must be contained within the TID region
    # (use all raw same-orientation breakpoints to bound the region)
    ext_pos = ([d[1] for d in mm_raw + pp_raw] + [d[4] for d in mm_raw + pp_raw])
    ext_lpos, ext_rpos = min(ext_pos), max(ext_pos)
    n_circ = sum(1 for cycle in cycleList if 0 not in cycle)
    bad = []
    for cycle in cycleList:
        if 0 in cycle:
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
