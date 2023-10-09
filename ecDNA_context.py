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
            bed_regions.append((chrom, int(start), int(end)))
    return bed_regions

def filter_graph_with_bed(graph_file_path, bed_regions):
    filtered_sequences = []
    filtered_edges = []
    sequence_endpoints = set()

    seq_count = 0
    
    with open(graph_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("sequence"):
                seq_count += 1
                _, start, end, cn, *_ = line.split("\t")
                start_chr, start_pos = start.split(":")[0], int(start.split(":")[1].split("-")[0])
                end_chr, end_pos = end.split(":")[0], int(end.split(":")[1].split("+")[0])
                
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
                
                if (start_chr, start_pos_int) in sequence_endpoints and (end_chr, end_pos_int) in sequence_endpoints:
                    filtered_edges.append({"type": type, "start_chr": start_chr, "end_chr": end_chr, "start": start_pos, "end": end_pos, "cn": float(cn)})
    return filtered_sequences, filtered_edges
 



def spanningChrs(filtered_sequences):
    #given filtered sequences, get the chromosomes involved in the amplicon
    chroms = set()
    for seq in filtered_sequences:
        start_chr, end_chr = seq['start_chr'], seq['end_chr']
        chroms.add(start_chr)
        chroms.add(end_chr)
    return list(chroms)

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
    
    print(f"cycle frac rows:\n{res}")
    return pd.DataFrame(res, columns=['Cycle','Size','CN','Circular','CN Fraction']).sort_values(by='CN Fraction', ascending=False)

def deduce_states(all_cns):
    #all_round = sorted(list(map(lambda x: round(x), all_cns)))
    all_round = sorted(all_cns)
    gaps = [all_round[i+1] - all_round[i] for i in range(len(all_round)-1)]
    if not gaps:
        return {all_round[0]: 1}
    threshold = 1
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


def count_transitions_del(filtered_sequences, filtered_edges, fuse_cutoff = 2000):
    discords = []
    for edge in filtered_edges:
        if edge['type'] == 'discordant':

            s_sign = edge['start'][-1]
            e_sign = edge['end'][-1]
            schr = edge['start_chr']
            echr = edge['end_chr']
            spos = pos_str_int(edge['start'])
            epos = pos_str_int(edge['end'])
            
            if e_sign == '-' and s_sign == '+' and schr == echr and spos > epos:
                #treat small deletions (less than fuse_cutoff) as artefacts and fuse together the sequences before and after
                if spos - epos < fuse_cutoff:
                    discords.append((epos, spos))
            elif e_sign == '+' and s_sign == '-' and schr == echr and spos > epos:
                if spos - epos < fuse_cutoff:
                    discords.append((epos + 1, spos - 1))
                    
    seqs_res = []
    fused = False
    for cur_seq in filtered_sequences:
        s, e = int(cur_seq['start']), int(cur_seq['end'])
        if (s,e) in discords:
            #fuse together previous and next sequences
            prev, later = cur_seq['ID'] - 1, cur_seq['ID'] + 1
            prev_seq = [seq for seq in filtered_sequences if seq['ID'] == prev][0]
            later_seq = [seq for seq in filtered_sequences if seq['ID'] == later][0]

            # fuse if prev and later have similar cn
            if abs(prev_seq['cn'] - later_seq['cn']) > 1 or prev_seq['end_chr'] != later_seq['start_chr']:
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
        
        if prev['end_chr'] == cur['start_chr'] and cur['start'] <= prev['end'] + 50:
            if abs(prev['cn'] - cur['cn']) > 1:
                trans += 1
       
    print(f"trans : {trans}")
    return trans, [seq['cn'] for seq in seqs_res] #return total transitions and all CNs      

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
 
def fbCount(filtered_edges, fb_cutoff = 25000):
    # count the number of foldback edges
    fb_count = 0
    for edge in filtered_edges:
        if edge['type'] == 'discordant':
            s_sign = edge['start'][-1]
            e_sign = edge['end'][-1]
            schr = edge['start_chr']
            echr = edge['end_chr']
            spos = pos_str_int(edge['start'])
            epos = pos_str_int(edge['end'])
            if e_sign == s_sign and schr == echr and abs(spos - epos) < fb_cutoff:
                fb_count += 1
    return fb_count

def detectTwoFB(filtered_sequences, filtered_edges):
    amp_start_chr, amp_start_pos = filtered_sequences[0]['start_chr'], filtered_sequences[0]['start']
    amp_end_chr, amp_end_pos = filtered_sequences[-1]['end_chr'], filtered_sequences[-1]['end']
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
                spos = pos_str_int(edge['start'])
                epos = pos_str_int(edge['end'])
                if e_sign == s_sign and schr == echr and schr == amp_start_chr:
                    if abs(spos - amp_start_pos) <= 50000:
                        start_found = True
                    elif abs(spos - amp_end_pos) <= 50000:
                        end_found = True
    return start_found and end_found

def ecDNAMetrics(filtered_cycles, filtered_sequences, filtered_edges, fuse_cutoff = 2000):
    cycle_frac = 0
    t_n_ratio = 0

    frac_table = cycleFractionTable(filtered_cycles, filtered_sequences)
    
    circular_frac = frac_table[frac_table['Circular']]
    sum_3_frac = 0.0
    if circular_frac.shape[0] > 0:
        cycle_frac = circular_frac['CN Fraction'].iloc[0]
        sum_3_frac = sum(circular_frac[:3]['CN Fraction'].values.tolist())
    
    n_trans, all_cns = count_transitions_del(filtered_sequences, filtered_edges, fuse_cutoff=fuse_cutoff)
    n_cn_states = len(deduce_states(all_cns))
    
    t_n_ratio = n_trans/n_cn_states
    
    n_chrs = len(spanningChrs(filtered_sequences))
    
    #cns_dict = trans_cn(gfile)
    cns_dict = deduce_states(all_cns)
    fb_count = fbCount(filtered_edges)

    two_fb = detectTwoFB(filtered_sequences, filtered_edges)
    
    return (cycle_frac, t_n_ratio, n_chrs, sum_3_frac, fb_count, len(filtered_edges), two_fb), n_trans, n_cn_states, cns_dict

def ecDNAContext(metrics, n_cn, t_n_cutoff = 4, cycle_cutoff = 0.15):
    cycle_frac, t_n_ratio, n_chrs, sum_3_frac, fb_count, tot_edges, two_fb = metrics
    if fb_count >= 4 and fb_count/tot_edges >= 0.1 and n_cn > 1:
        return "BFB"
    
    elif t_n_ratio < t_n_cutoff and two_fb:
        return "Two-foldback"
    
    elif cycle_frac >= cycle_cutoff and (n_cn == 1):
        return "Simple circular"
    
    elif cycle_frac >= cycle_cutoff and t_n_ratio >= t_n_cutoff and n_cn > 1:
        if n_chrs > 1:
            return "Both, multi-chromosomal"
        else:
            return "Both, single-chromosomal"
    
    elif t_n_ratio >= t_n_cutoff and n_cn > 1:
        if n_chrs > 1:
            return "Heavily Rearranged, Multiple chromosomes"
        else:
            return "Heavily Rearranged"
    else:
        return "Unknown"
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the context behind ecDNA formation")
    parser.add_argument("-g", "--graph", help="Path to amplicon graph file", required=True)
    parser.add_argument("-c", "--cycles", help="Path to cycles file", required=True)
    parser.add_argument("-f", "--fuse", help="Fuse cutoff for deletions", default=2000)
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
    filtered_sequences, filtered_edges = filter_graph_with_bed(graph_file, bed_regions)
    _, converted_cycles = make_new_cycle(graph_file, cycles_file)
    filtered_cycles = filter_cycles_with_edges(converted_cycles, filtered_sequences)

    metrics, n_trans, n_cn, cn_dict = ecDNAMetrics(filtered_cycles, filtered_sequences, filtered_edges, FUSE_CUTOFF)
    context = ecDNAContext(metrics, n_cn, cycle_cutoff=CYCLE_FRAC_CUTOFF, t_n_cutoff=TN_RATIO_CUTOFF)

    # output all the information
    print(f"Metrics: {metrics}")
    print(f"n_trans: {n_trans}")
    print(f"n_cn: {n_cn}")
    print(f"cn_dict: {cn_dict}")
    print(f"Context: {context}")

