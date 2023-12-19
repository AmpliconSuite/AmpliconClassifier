#!/usr/bin/env python3

import os
import argparse

#Author: Siavash Dehkordi, modified by Jens Luebeck and Bhargavi Dameracharla

#converts cycles file to be numbered by BPG segment.

def extract_seg(path):
    line_number = 0
    seg_dict = {}
    name = ''
    segments_list = [{
                    'segment_number': 0,
                    'start_chr': None,
                    'start_index': None,
                    'end_chr': None,
                    'end_index': None,
                    'cn': None
                }]
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line_number == 0:
                name = ''
                # continue
                # name = line.split(' ')[3]
            elif line.startswith('sequence'):
                line = line.split()
                segment_number = line_number
                line_s = line[1].split(':')
                start_chr = line_s[0]
                start_index = line_s[1][:-1]
                line_e = line[2].split(':')
                end_chr = line_e[0]
                end_index = line_e[1][:-1]
                cn = line[3]
                seg_dict[start_index] = {'segment_number': segment_number, 'start_chr': start_chr,
                                         'start_index': start_index,
                                         'end_cht': end_chr, 'end_index': end_index}
                seg_dict[end_index] = {'segment_number': segment_number, 'start_chr': start_chr,
                                       'start_index': start_index,
                                       'end_cht': end_chr, 'end_index': end_index}
                segments_list.append({
                    'segment_number': segment_number,
                    'start_chr': start_chr,
                    'start_index': start_index,
                    'end_chr': end_chr,
                    'end_index': end_index,
                    'cn': cn
                })
            line_number += 1
    return name, seg_dict, segments_list


def extract_fuse_seg(path):
    fuse_seg_dict = {}
    with open(path) as f:
        for line in f:
            if line.startswith('Segment'):
                line = line.strip()
                line = line.split()
                fuse_seg_dict[line[1]] = {'start_chr': line[2], 'start_index': line[3],
                                          'end_cht': line[2], 'end_index': line[4]}
    return fuse_seg_dict


def make_new_cycle(graph_path, cycle_path):
    _, seg_dict, segments_list = extract_seg(graph_path)
    fuse_seg_dict = extract_fuse_seg(cycle_path)
    segment_dict = {'0': {'segment_chr': '0', 'segment_start': '0',
                          'segment_end': '0', 'length': '0'}}

    converted_cycles = []
    with open(cycle_path) as f:
        for line in f:
            if line.startswith('Segment'):
                line = line.strip().split()
                segment_number = line[1]
                segment_chr = line[2]
                segment_start = line[3]
                segment_end = line[4]
                segment_length = int(segment_end) - int(segment_start)
                segment_dict[segment_number] = {'segment_chr': segment_chr, 'segment_start': segment_start,
                                                'segment_end': segment_end, 'length': segment_length}
            elif line.startswith('Cycle'):
                cycle_info = {}
                continues = 1
                end_loc = 0
                length = 0
                fields = line.split(';')
                splitfields = [(x.rsplit('=')[0], x.rsplit('=')[1]) for x in fields]
                for x in splitfields:
                    if x[0] != "Segments":
                        cycle_info[x[0]] = x[1]
                    else:
                        segseq = x

                a = segseq[0] + '='
                new_segs = []
                circle = True
                segseq = segseq[1].strip().split(',')
                for l in segseq:
                    number = l[:-1]
                    length = length + int(segment_dict[number]['length'])
                    sign = l[-1]
                    if number == '0':
                        circle = False
                        new_segs.append(l)
                    else:
                        fuse_seg = fuse_seg_dict[number]
                        start_seg = int(fuse_seg['start_index'])
                        end_seg = int(fuse_seg['end_index'])
                        start_number = int(seg_dict[str(start_seg)]['segment_number'])
                        end_number = int(seg_dict[str(end_seg)]['segment_number'])
                        if sign == '+':
                            for i in range(start_number, end_number + 1):
                                new_segs.append(str(i) + sign)
                            if end_loc == 0:
                                end_loc = int(segment_dict[number]['segment_end'])
                                chrom = segment_dict[number]['segment_chr']
                            else:
                                if (chrom != segment_dict[number]['segment_chr'] or abs(
                                        int(segment_dict[number]['segment_end']) - end_loc) > 10) and number != 0:
                                    continues +=1
                                end_loc = int(segment_dict[number]['segment_end'])
                                chrom = segment_dict[number]['segment_chr']
                        elif sign == '-':
                            for i in range(end_number, start_number - 1, -1):
                                new_segs.append(str(i) + sign)
                           
                            if end_loc == 0:
                                end_loc = int(segment_dict[number]['segment_start'])
                                chrom = segment_dict[number]['segment_chr']
                            else:
                                if (chrom != segment_dict[number]['segment_chr'] or abs(
                                        int(segment_dict[number]['segment_end']) - end_loc) > 10) and number != 0:
                                    continues +=1
                                end_loc = int(segment_dict[number]['segment_end'])
                                chrom = segment_dict[number]['segment_chr']
                
                cycle_info['Length'] = str(int(length / 1000)) + 'Kbp'
                cycle_info['Circular'] = circle
                cycle_info['Segments'] = ','.join(new_segs)
                if continues > 2:
                    cycle_info['TRIVIAL'] = False
                else:
                    cycle_info['TRIVIAL'] = True

                converted_cycles.append(cycle_info)
    return segments_list, converted_cycles
