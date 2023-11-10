import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from  collections  import  defaultdict
from convert_cycles_file import make_new_cycle
import argparse

def read_bed_file(bed_file_path):
    bed_regions = []
    with open(bed_file_path, 'r') as rf:
        for line in rf:
            chrom, start, end = line.strip().split()[:3]
            bed_regions.append((chrom, int(start)-10000, int(end)+10000))
    return bed_regions

def regions_apart(bed_regions):
    #check if there are any pair of segments are far apart (either different chromosomes  or more than 1 MB apart)
    for i in range(len(bed_regions)-1):
        for j in range(i+1, len(bed_regions)):
            if bed_regions[i][0] != bed_regions[j][0] or abs(bed_regions[j][1] - bed_regions[i][2]) > 1000000:
                return True
    return False

def filter_graph_with_bed(graph_file_path, bed_regions):
    filtered_sequences = []
    filtered_edges = []
    sequence_endpoints = set()

    all_sequences = []
    all_edges = []

    seq_count = 0
    
    with open(graph_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("sequence"):
                seq_count += 1
                _, start, end, cn, *_ = line.split("\t")
                start_chr, start_pos = start.split(":")[0], int(start.split(":")[1].split("-")[0])
                end_chr, end_pos = end.split(":")[0], int(end.split(":")[1].split("+")[0])

                all_sequences.append({"ID": seq_count,"start_chr": start_chr, "end_chr": end_chr, "start": start_pos, "end": end_pos, "cn": float(cn)})

                for chrom, bed_start, bed_end in bed_regions:
                    if chrom == start_chr and bed_start <= end_pos and bed_end >= start_pos:
                        filtered_sequences.append({"ID": seq_count,"start_chr": start_chr, "end_chr": end_chr, "start": start_pos, "end": end_pos, "cn": float(cn)})
                        sequence_endpoints.add((start_chr, start_pos))
                        sequence_endpoints.add((end_chr, end_pos))
                        break

                        
            elif line.startswith("source") or line.startswith("concordant") or line.startswith("discordant"):
                type, edge_info, cn, *_ = line.split("\t")
                start_info, end_info = edge_info.split("->")
                start_chr, start_pos = start_info.split(":")
                start_pos_int = int(start_pos[:-1])
                end_chr, end_pos = end_info.split(":")
                end_pos_int = int(end_pos[:-1])

                all_edges.append({"type": type, "start_chr": start_chr, "end_chr": end_chr, "start": start_pos, "end": end_pos, "cn": float(cn)})
                if (start_chr, start_pos_int) in sequence_endpoints and (end_chr, end_pos_int) in sequence_endpoints:
                    filtered_edges.append({"type": type, "start_chr": start_chr, "end_chr": end_chr, "start": start_pos, "end": end_pos, "cn": float(cn)})

    return filtered_sequences, filtered_edges, all_sequences, all_edges



def spanningChrs(filtered_sequences):
    #given filtered sequences, get the chromosomes involved in the amplicon
    chroms = set()
    chrom_contrib = defaultdict(int)
    for seq in filtered_sequences:
        start_chr, end_chr, start, end = seq['start_chr'], seq['end_chr'], seq['start'], seq['end']
        chroms.add(start_chr)
        chrom_contrib[start_chr] += abs(int(start) - int(end))
    return list(chroms), chrom_contrib

def seg_str_int(segs):
    res = []
    for seg in segs:
        i = int(seg[:len(seg)-1])
        if seg[-1] == '-': i *= -1
        res.append(i)
    return res

def pos_str_int(pos):
    #get -i or +i int from str i- or i+
    res = int(pos[:len(pos)-1])
    if pos[-1] == '-':
        res *= -1
    return res

def cycleOverlapFraction(edge, filtered_sequences):
    #Get the fraction of copy numbee accounted by sequences that overlap with the cycle (single everted edge)
    start, end = abs(pos_str_int(edge['end'])), abs(pos_str_int(edge['start']))
    # end, start flipped above to maintain increasing order
    total_cn = 0
    overlap_cn = 0
    for seq in filtered_sequences:
        s, e = int(seq['start']), int(seq['end'])
        cn = float(seq['cn'])
        total_cn += cn
        #check overlap - has to lie fully within cycle
        if s >= start and e <= end:
            overlap_cn += cn
    return overlap_cn/total_cn

def cycleFractionTable(filtered_cycles, filtered_sequences):
    res = []
    total_cn = 0 #get total weighted copy number for the amplicon
    # segments are not overlapping
    for cur_seg in filtered_sequences:
        chrom, start, end, seg_cn = cur_seg['start_chr'], int(cur_seg['start']), int(cur_seg['end']), float(cur_seg['cn'])
        total_cn += seg_cn*(abs(end - start))

    for c in filtered_cycles:
        cycID = c['Cycle']
        cycle_cn, cycle_size = float(c['Copy_count']), c['Length']
        circular = c['Circular']
        cycle_size_num = int(cycle_size.split("Kbp")[0])*1000
        cycle_frac = cycle_size_num*cycle_cn/total_cn
        row = [cycID, cycle_size,cycle_cn, circular, cycle_frac]
        res.append(row)
    
    return pd.DataFrame(res, columns=['Cycle','Size','CN','Circular','CN Fraction']).sort_values(by='CN Fraction', ascending=False)

def deduce_states(all_cns):
    #all_round = sorted(list(map(lambda x: round(x), all_cns)))
    all_round = sorted(all_cns)
    gaps = [all_round[i+1] - all_round[i] for i in range(len(all_round)-1)]
    if not gaps:
        return {all_round[0]: 1}
    threshold = 2
    cur_cns = []
    count_dict = {}
    for i in range(len(all_round)):
        if len(cur_cns) > 0:
            if abs(np.mean(cur_cns) - all_round[i]) >= threshold:
                count_dict[round(sum(cur_cns)/len(cur_cns))] = len(cur_cns)
                cur_cns = []
        cur_cns.append(all_round[i])
    count_dict[round(sum(cur_cns)/len(cur_cns))] = len(cur_cns)
    return count_dict


def count_transitions_del(filtered_sequences, filtered_edges, fuse_cutoff = 5000):
    discords = []
    n_dels = 0
    for edge in filtered_edges:
        if edge['type'] == 'discordant':

            s_sign = edge['start'][-1]
            e_sign = edge['end'][-1]
            schr = edge['start_chr']
            echr = edge['end_chr']
            spos = abs(pos_str_int(edge['start']))
            epos = abs(pos_str_int(edge['end']))
            
            if e_sign == '-' and s_sign == '+' and schr == echr and spos > epos:
                #treat small deletions (less than fuse_cutoff) as artefacts and fuse together the sequences before and after
                if spos - epos < fuse_cutoff:
                    discords.append((epos, spos))
                else:
                    n_dels += 1
            elif e_sign == '+' and s_sign == '-' and schr == echr and spos > epos:
                if spos - epos < fuse_cutoff:
                    discords.append((epos + 1, spos - 1))
                else:
                    n_dels += 1
                    
    seqs_res = []
    fused = False
    for cur_seq in filtered_sequences:
        s, e = int(cur_seq['start']), int(cur_seq['end'])
        if (s,e) in discords:
            #fuse together previous and next sequences
            prev, later = cur_seq['ID'] - 1, cur_seq['ID'] + 1
            prev_seq = [seq for seq in filtered_sequences if seq['ID'] == prev]

            later_seq = [seq for seq in filtered_sequences if seq['ID'] == later]
            if not prev_seq or not later_seq:
                continue
            else:
                prev_seq = prev_seq[0]
                later_seq = later_seq[0]

            # fuse if prev and later have similar cn
            if abs(prev_seq['cn'] - later_seq['cn']) > 2 or prev_seq['end_chr'] != later_seq['start_chr']:
                seqs_res.append(cur_seq)
                continue
            else:
                new_seq = {'ID':prev_seq['ID'],'start_chr':prev_seq['start_chr'],'end_chr':later_seq['end_chr'],'start':prev_seq['start'],'end':later_seq['end'],'cn':prev_seq['cn'],'seg_size':later_seq['end'] - prev_seq['start']}
                seqs_res.append(new_seq)
                fused = True
        else:
            #check if current sequence was already added
            if fused:
                fused = False
                continue
            else:
                seqs_res.append(cur_seq)
                
    #now count transitions
    trans = 0
    for i in range(1, len(seqs_res)):
        prev = seqs_res[i-1]
        cur = seqs_res[i]
        
        #if prev['end_chr'] == cur['start_chr'] and cur['start'] <= prev['end'] + 50:
        if (abs(prev['cn'] - cur['cn']) > 2 and (abs(prev['end'] - prev['start'] > 5000) and abs(cur['end'] - cur['start']) > 5000 )) or cur['start'] - prev['end'] > 5000:
            trans += 1
       
    return trans, [seq['cn'] for seq in seqs_res if abs(seq['end'] - seq['start']) > 5000]

def filter_cycles_with_edges(cycles, filtered_sequences):
    # keep the cycles where all the segments in the cycle are in filtered_sequences
    filtered_cycles = []
    allowed_segs = set()
    allowed_segs.add(0)
    for seq in filtered_sequences:
        allowed_segs.add(int(seq['ID']))
    for cycle in cycles:
        cur_segs = cycle['Segments'].split(',')
        cur_segs = [seg[:-1] for seg in cur_segs]
        cur_segs = [int(seg) for seg in cur_segs]
        if set(cur_segs).issubset(allowed_segs):
            filtered_cycles.append(cycle)
    return filtered_cycles
 
def fbCount(filtered_edges, fb_cutoff = 50000):
    # count the number of foldback edges
    fb_count = 0
    for edge in filtered_edges:
        if edge['type'] == 'discordant':
            s_sign = edge['start'][-1]
            e_sign = edge['end'][-1]
            schr = edge['start_chr']
            echr = edge['end_chr']
            spos = abs(pos_str_int(edge['start']))
            epos = abs(pos_str_int(edge['end']))
            if e_sign == s_sign and schr == echr and abs(spos - epos) < fb_cutoff:
                fb_count += 1
    return fb_count

def detectTwoFB(filtered_sequences, filtered_edges, fb_cutoff = 50000):
    amp_start_chr, amp_start_pos = filtered_sequences[0]['start_chr'], filtered_sequences[0]['start']
    amp_start_end = filtered_sequences[0]['end']
    amp_end_chr, amp_end_pos = filtered_sequences[-1]['end_chr'], filtered_sequences[-1]['end']
    amp_end_start = filtered_sequences[-1]['start']
    start_found = False
    end_found = False
    if amp_start_chr == amp_end_chr:
        # look for foldback edges closer to these end points
        for edge in filtered_edges:
            if edge['type'] == 'discordant':
                s_sign = edge['start'][-1]
                e_sign = edge['end'][-1]
                schr = edge['start_chr']
                echr = edge['end_chr']
                spos = abs(pos_str_int(edge['start']))
                epos = abs(pos_str_int(edge['end']))
                if abs(epos - spos) <= fb_cutoff:
                    if e_sign == s_sign and schr == echr and schr == amp_start_chr:
                        if e_sign == '-':
                            if abs(spos - amp_start_pos) <= 50000 or abs(epos - amp_start_pos) <= 50000 or abs(spos - amp_start_end) <= 50000 or abs(epos - amp_start_end) <= 50000:
                                start_found = True
                        else:
                            if abs(spos - amp_end_pos) <= 50000 or abs(epos - amp_end_pos) <= 50000 or abs(spos - amp_end_start) <= 50000 or abs(epos - amp_end_start) <= 50000:
                                end_found = True
    return start_found and end_found

def crossEdges(filtered_edges):
    # Among the discordant edges, how many cross each other?
    # Two edges cross when one of the breakpoints of one edge is between the breakpoints of the other edge
    # Only look at discordant edges
    cross_count = 0
    for i in range(len(filtered_edges)-1):
        edge1 = filtered_edges[i]
        if edge1['type'] != 'discordant':
            continue
        s1_chr = edge1['start_chr']
        e1_chr = edge1['end_chr']
        s1, e1 = abs(pos_str_int(edge1['start'])), abs(pos_str_int(edge1['end']))
        s1, e1 = min(s1, e1), max(s1, e1)

        if s1_chr != e1_chr:
            cross_count += 1

        for j in range(i+1, len(filtered_edges)):
            edge2 = filtered_edges[j]
            if edge2['type'] != 'discordant':
                continue
            s2_chr = edge2['start_chr']
            e2_chr = edge2['end_chr']
            s2, e2 = abs(pos_str_int(edge2['start'])), abs(pos_str_int(edge2['end']))
            s2, e2 = min(s2, e2), max(s2, e2)

            if j == len(filtered_edges) - 1:
                if s2_chr != e2_chr:
                    cross_count += 1

            if s1_chr == e1_chr:
                if s2_chr == s1_chr:
                    if s1 < s2 < e1 and e2 > e1:
                        cross_count += 1
                        continue
                elif e2_chr == s1_chr:
                    if s1 < e2 < e1 and s2 < s1:
                        cross_count += 1
                        continue
            if s2_chr == e2_chr:
                if s1_chr == s2_chr:
                    if s2 < s1 < e2 and e1 > e2:
                        cross_count += 1
                        continue
                elif e1_chr == s2_chr:
                    if s2 < e1 < e2 and s1 < s2:
                        cross_count += 1
                        continue
            
    return cross_count

def small_contribution(chrom_contribs):
    # if there is only a small segment of a chromosome, we allow it to be ignored when checking for episome
    chrom_lengths = [chrom_contribs[c] for c in chrom_contribs if chrom_contribs[c] > 30000]
    if len(chrom_lengths) >= 2:
        return False
    return True
            
def ecDNAMetrics(filtered_cycles, filtered_sequences, filtered_edges, sequences, edges, multiple_large, fuse_cutoff = 5000):
    
    two_fb = detectTwoFB(filtered_sequences, filtered_edges)
    n_chrs, chrom_contribs = spanningChrs(filtered_sequences)
    n_chrs = len(n_chrs)
    small_contrib = small_contribution(chrom_contribs)

    cycle_frac = 0
    t_n_ratio = 0
    frac_table = cycleFractionTable(filtered_cycles, filtered_sequences)
    circular_frac = frac_table[frac_table['Circular']]
    
    use_everted_metrics = False
    if circular_frac.shape[0] > 0:
        cycle_frac = circular_frac['CN Fraction'].iloc[0]
        largest_cycle_id = circular_frac['Cycle'].iloc[0]
        largest_cycle = [c for c in filtered_cycles if c['Cycle'] == largest_cycle_id][0]

        everted, outside_seqs, outside_edges = simpleEverted(largest_cycle, filtered_sequences, filtered_edges)
        if everted and not multiple_large:
            n_trans, all_cns = count_transitions_del(outside_seqs, outside_edges, fuse_cutoff=fuse_cutoff)
            n_cn_states = len(deduce_states(all_cns))
            t_n_ratio = n_trans/n_cn_states
            fb_count = fbCount(outside_edges)
            cross_edges = crossEdges(outside_edges)

            # check if the rest of the amplicon is simple enough to be considered episome even if there is rearrangement inside the cycle.
            if (cross_edges <= 1) and n_cn_states <= 3 and fb_count <= 3 and (n_chrs == 1 or small_contrib):
                use_everted_metrics = True

    if not use_everted_metrics:
        n_trans, all_cns = count_transitions_del(filtered_sequences, filtered_edges, fuse_cutoff=fuse_cutoff)
        n_cn_states = len(deduce_states(all_cns))
        
        t_n_ratio = n_trans/n_cn_states
        
        cns_dict = deduce_states(all_cns)
        fb_count = fbCount(filtered_edges)
        cross_edges = crossEdges(filtered_edges)

        # to check if the rest of the amplicon is complex, we need to look at all sequences, even those outside ecDNA region
        sequences2 = sequences
        edges2 = edges
    
    else:
        # Look at all sequences in the amplicon to understand if ecDNA occurs with a complex background but if we know that it is a simple cycle with internal rearrangement, remove the rearrangements inside the cycle.
        _, sequences2, edges2 = simpleEverted(largest_cycle, sequences, edges)

    
    n_trans2, all_cns2 = count_transitions_del(sequences2, edges2, fuse_cutoff=fuse_cutoff)
    n_cn_states2 = len(deduce_states(all_cns2))
    t_n_ratio2 = n_trans2/n_cn_states2
    cross_edges2 = crossEdges(edges2)

    #check if only one discordant edge and get overlap fraction. There can be other concordant edges
    disc_edges = [edge for edge in filtered_edges if edge['type'] == 'discordant']
    if len(disc_edges) == 1:
        edge = disc_edges[0]
        overlap_frac = cycleOverlapFraction(edge, filtered_sequences)
    else:
        overlap_frac = 0
    
    metrics = {}
    metrics['cycle_frac'] = cycle_frac
    metrics['t_n_ratio'] = t_n_ratio
    metrics['n_chrs'] = n_chrs
    metrics['n_cn'] = n_cn_states
    metrics['n_trans'] = n_trans
    metrics['fb_count'] = fb_count
    metrics['total_edges'] = len(filtered_edges)
    metrics['two_fb'] = two_fb
    metrics['cross_edges'] = cross_edges
    metrics['small_contrib'] = small_contrib
    metrics['overlap_frac'] = overlap_frac
    metrics['use_everted_metrics'] = use_everted_metrics

    metrics['n_trans_amplicon'] = n_trans2
    metrics['n_cn_amplicon'] = n_cn_states2
    metrics['t_n_ratio_amplicon'] = t_n_ratio2
    metrics['cross_edges_amplicon'] = cross_edges2

    return metrics

# if the biggest circle is formed due to an everted edge, we don't need to consider the rerarrangements/copy number states within. It should be called simple circular in this case. 
def simpleEverted(largest_cycle, filtered_sequences, filtered_edges):
    cycle_segs = largest_cycle['Segments'].split(",")
    cycle_start_seg, cycle_end_seg = cycle_segs[0], cycle_segs[-1]
    #get the start of start seq and end of end seq
    start_seq = None
    end_seq = None
    for seq in filtered_sequences:
        if seq['ID'] == abs(pos_str_int(cycle_start_seg)):
            start_seq = seq
        if seq['ID'] == abs(pos_str_int(cycle_end_seg)):
            end_seq = seq
    edge_start = end_seq['end']
    edge_end = start_seq['start']

    #look for everted edge that is edge_start+ -> edge_end-
    everted = None
    edge_found = False
    if start_seq['start_chr'] == end_seq['end_chr']:
        for e in filtered_edges:
            if e['type'] == 'discordant' and e['start_chr'] == start_seq['start_chr'] and e['end_chr'] == e['start_chr'] and e['start'] == str(edge_start) + '+' and e['end'] == str(edge_end) + '-':
                everted = e
                edge_found = True
    
    if edge_found:
        #remove all sequences, edges that fall within these two positions and get all the metrics
        sequences_outside = []
        edges_outside = []
        for seq in filtered_sequences:
            if seq['start'] <= edge_end or seq['end'] >= edge_start or seq['start_chr'] != everted['start_chr']:
                sequences_outside.append(seq)

        for e in filtered_edges:
            if e['start_chr'] != everted['start_chr'] or e['end_chr'] != everted['end_chr'] or abs(pos_str_int(e['start'])) <= edge_end or abs(pos_str_int(e['end'])) >= edge_start:
                edges_outside.append(e)
        return everted, sequences_outside, edges_outside
    
    return edge_found, filtered_sequences, filtered_edges



def ecDNAContext(metrics, t_n_cutoff = 4, cycle_cutoff = 0.15):
    n_cn = metrics['n_cn']
    t_n_ratio = metrics['t_n_ratio']
    n_cross_edges = metrics['cross_edges']
    n_chrs = metrics['n_chrs']
    cycle_frac = metrics['cycle_frac']
    small_contrib = metrics['small_contrib']

    if (metrics['fb_count'] >= 4 and (metrics['fb_count']/metrics['total_edges'] >= 0.1)) and n_cn > 2 and n_chrs == 1:
        return "BFB-like"
    
    elif ((t_n_ratio < t_n_cutoff) or n_cross_edges < 4) and metrics['two_fb'] and n_cn <= 6:
        return "Two-foldback"

    elif ((t_n_ratio >= t_n_cutoff or n_cross_edges >= 4) and n_cn > 1) or (n_cn > 6):
        if n_chrs > 1 and not small_contrib:
            return "Heavily rearranged, multi-chromosomal"
        else:
            return "Heavily rearranged"
    
    elif cycle_frac >= cycle_cutoff and ((t_n_ratio < t_n_cutoff) or n_cross_edges < 4) and n_cn <= 6 and (n_chrs == 1 or small_contrib):
        if metrics['t_n_ratio_amplicon'] >= t_n_cutoff or metrics['cross_edges_amplicon'] >= 4 or metrics['n_cn_amplicon'] > 6:
            return "Simple circular, complex background"
        return "Simple circular"
    
    elif metrics['overlap_frac'] >= 0.3:
        if metrics['t_n_ratio_amplicon'] >= t_n_cutoff or metrics['cross_edges_amplicon'] >= 4 or metrics['n_cn_amplicon'] > 6:
            return "Simple circular, complex background"
        return "Simple circular"
    else:
        return "Unknown"
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the context behind ecDNA formation")
    parser.add_argument("-g", "--graph", help="Path to amplicon graph file", required=True)
    parser.add_argument("-c", "--cycles", help="Path to cycles file", required=True)
    parser.add_argument("-f", "--fuse", help="Fuse cutoff for deletions", default=5000)
    parser.add_argument("-t", "--t_n", help="Transition to CN ratio cutoff", default=4)
    parser.add_argument("-y", "--cycle_cutoff", help="Cycle fraction cutoff", default=0.15)
    parser.add_argument("-b", "--bed", help="Bed file of ecDNA regions", required=False)

    args = parser.parse_args()
    graph_file = args.graph
    cycles_file = args.cycles
    FUSE_CUTOFF = args.fuse
    TN_RATIO_CUTOFF = args.t_n
    CYCLE_FRAC_CUTOFF = args.cycle_cutoff
    bed_file = args.bed

    bed_regions = read_bed_file(bed_file)
    multiple_large = regions_apart(bed_regions)
    filtered_sequences, filtered_edges, sequences, edges = filter_graph_with_bed(graph_file, bed_regions)
    _, converted_cycles = make_new_cycle(graph_file, cycles_file)
    filtered_cycles = filter_cycles_with_edges(converted_cycles, filtered_sequences)

    metrics= ecDNAMetrics(filtered_cycles, filtered_sequences, filtered_edges, sequences, edges, multiple_large, FUSE_CUTOFF)
    context = ecDNAContext(metrics, cycle_cutoff=CYCLE_FRAC_CUTOFF, t_n_cutoff=TN_RATIO_CUTOFF)

    # output all the information
    print(f"Metrics: {metrics}")
    print(f"Context: {context}")

