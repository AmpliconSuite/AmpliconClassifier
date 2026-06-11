#!/usr/bin/env python3
"""
Find non-SV-spanned CN drop regions consistent across amplicons.
Outputs an HTML file with clickable AA image links for representative examples.

Usage: python find_patch_candidates.py input_file.input [out_html] [cn_ratio] [min_frac] [min_n]
"""
import sys, os, collections, statistics

def parse_graph(gf):
    segs, disc = [], set()
    in_seq = False
    with open(gf) as f:
        for line in f:
            t = line.strip(); p = t.split('\t')
            if t.startswith('SequenceEdge:'):   in_seq = True;  continue
            if t.startswith('BreakpointEdge:'): in_seq = False; continue
            if in_seq and len(p) >= 4:
                try:
                    c = p[1].split(':')[0]
                    pl = int(p[1].split(':')[1][:-1]); pr = int(p[2].split(':')[1][:-1])
                    segs.append((c, min(pl, pr), max(pl, pr), float(p[3])))
                except: pass
            elif not in_seq and p[0] == 'discordant' and len(p) >= 2:
                try:
                    for v in p[1].split('->'):
                        disc.add((v.split(':')[0], int(v.split(':')[1][:-1])))
                except: pass
    return sorted(segs), disc

def main():
    infile   = sys.argv[1]
    outhtml  = sys.argv[2] if len(sys.argv) > 2 else infile.replace('.input', '_patch_candidates.html')
    ratio    = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
    min_frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
    min_n    = int(sys.argv[5])   if len(sys.argv) > 5 else 20
    min_flank_cn = 4.5
    bin_sz   = 1000

    drop_c    = collections.Counter()
    cov_c     = collections.Counter()
    drop_data = collections.defaultdict(list)   # key -> [(gap_size, cn_ratio, sample, png)]
    N_REPS    = 50

    with open(infile) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 3: continue
            sample = p[0]
            gf = p[2]
            png = gf.replace('_graph.txt', '.png')
            try: segs, disc = parse_graph(gf)
            except: continue

            by_chrom = collections.defaultdict(list)
            for c, s, e, cn in segs:
                by_chrom[c].append((s, e, cn))

            amp_covered = set()
            amp_drops   = {}

            for chrom, cs in by_chrom.items():
                cs.sort()
                for i in range(1, len(cs) - 1):
                    s, e, cn = cs[i]
                    flank = (cs[i-1][2] + cs[i+1][2]) / 2
                    if flank < min_flank_cn: continue
                    key = (chrom, (s // bin_sz) * bin_sz, (e // bin_sz) * bin_sz)
                    amp_covered.add(key)
                    if (key not in amp_drops
                            and cn < ratio * flank
                            and (chrom, s) not in disc and (chrom, s - 1) not in disc
                            and (chrom, e) not in disc and (chrom, e + 1) not in disc):
                        amp_drops[key] = (e - s, cn / flank, sample, png)

            for key in amp_covered:
                cov_c[key] += 1
            for key, vals in amp_drops.items():
                drop_c[key] += 1
                drop_data[key].append(vals)

    # Build results, pick 5 reps closest to median CN ratio
    results = []
    for key in drop_c:
        nd, nc = drop_c[key], cov_c[key]
        if nd < min_n or nd / nc < min_frac: continue
        data = drop_data[key]
        sizes, ratios, samples, pngs = zip(*data)
        med_ratio = statistics.median(ratios)
        med_size  = int(statistics.median(sizes))
        reps = sorted(data, key=lambda x: abs(x[1] - med_ratio))[:N_REPS]
        results.append((key, nd, nc, med_size, med_ratio, reps))
    results.sort(key=lambda x: -x[1])

    # Write HTML
    with open(outhtml, 'w') as out:
        out.write('''<!DOCTYPE html><html><head><meta charset="utf-8">
<style>
  body { font-family: monospace; font-size: 13px; margin: 16px; }
  h2 { margin: 0 0 4px 0; font-size: 15px; }
  .stats { color: #555; margin-bottom: 8px; }
  .candidate { border: 1px solid #ccc; border-radius: 4px;
               padding: 10px; margin-bottom: 24px; }
  .gallery { display: flex; flex-wrap: wrap; gap: 8px; }
  .thumb { text-align: center; }
  .thumb img { max-height: 220px; width: auto; border: 1px solid #ddd;
               display: block; }
  .thumb span { font-size: 11px; color: #333; display: block;
                max-width: 200px; word-break: break-all; }
  .toc { margin-bottom: 24px; }
  .toc a { margin-right: 12px; color: #0055cc; }
</style></head><body>\n''')

        # Table of contents
        out.write('<div class="toc"><b>Candidates:</b><br>\n')
        for (chrom, bs, be), nd, nc, med_size, med_ratio, reps in results:
            anchor = f'{chrom}_{bs}_{be}'
            out.write(f'<a href="#{anchor}">{chrom}:{bs:,}–{be:,} '
                      f'(n={nd}, {nd/nc:.0%}, ratio={med_ratio:.3f})</a>\n')
        out.write('</div>\n')

        for (chrom, bs, be), nd, nc, med_size, med_ratio, reps in results:
            anchor = f'{chrom}_{bs}_{be}'
            out.write(f'<div class="candidate" id="{anchor}">\n')
            out.write(f'<h2>{chrom}:{bs:,}–{be:,}</h2>\n')
            out.write(f'<div class="stats">'
                      f'n={nd} &nbsp;|&nbsp; covered={nc} &nbsp;|&nbsp; '
                      f'fraction={nd/nc:.3f} &nbsp;|&nbsp; '
                      f'gap={med_size:,} bp &nbsp;|&nbsp; '
                      f'median CN ratio={med_ratio:.3f}</div>\n')
            out.write('<div class="gallery">\n')
            for gap, r, samp, png in reps:
                if os.path.exists(png):
                    out.write(f'<div class="thumb">'
                              f'<img src="file://{png}" alt="{samp}">'
                              f'<span>{samp}<br>ratio={r:.3f}</span>'
                              f'</div>\n')
            out.write('</div></div>\n')

        out.write('</body></html>\n')

    print(f'Wrote {len(results)} candidates to {outhtml}')

if __name__ == '__main__':
    main()
