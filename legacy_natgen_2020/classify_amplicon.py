
import sys, os, subprocess, copy
import hg19util as hg19

#Only include complete human chrs for classification
human_chrs = ["%d" % x for x in xrange(1,23)]
human_chrs.extend(['Y','X'])
human_chrs.extend(['chrY','chrX'])
human_chrs.extend(["chr%d" % x for x in xrange(1,23)])
decoys = set(['hs37d5','chrhs37d5'])

def classify_final_amplicon(amplicon, threshold_copy = 4, min_length = 10000):
  intervals = hg19.interval_list([i for i in amplicon['interval_map'].values()])
  intervals.sort()
  #Search for valid cycles
  cycles = amplicon['cycle_map'].keys()
  cycles = sorted(cycles, key=lambda x: amplicon['cycle_map'][x]['copy_count'], reverse=True)
  for c in amplicon['cycle_map']:
    #Cycles start with non-zero initial segment
    if amplicon['cycle_map'][c]['cycle'][0][0] != '0':
      #Check if any segments are from the decoy regions, if so, skip
      hits = len([s for s in amplicon['cycle_map'][c]['cycle'] if amplicon['segment_map'][s[0:-1]].chrom in decoys])
      if hits > 0:
        continue
      #Get all human segments that make up cycle          
      segments = hg19.interval_list([amplicon['segment_map'][s[0:-1]] for s in amplicon['cycle_map'][c]['cycle'] if amplicon['segment_map'][s[0:-1]].chrom in human_chrs])      
      segments.sort()
      #Find total length of cycle
      total_length = sum([s.end-s.start for s in segments])
      if total_length == 0:
        continue
      #Find total weight of cycle (length * copy number)
      wtotal_length = sum([(s.end-s.start)*s.info['AvAmp']*2 for s in segments])
      #Find average copy number of cycle
      copy = wtotal_length/total_length
      #As long as total length >= min_length and average copy >= min_copy, this is circular
      if total_length >= min_length and copy >= threshold_copy:
        return 'Circular'
  #Find total length of amplified region of amplicon, and if the total length is less than threshold, this is an invalid amplicon
  total_length = sum([s.end-s.start for s in amplicon['graph']['segments'] if s.info['Amplification']*2 >= threshold_copy and s.chrom in human_chrs])
  if total_length <= min_length:
    return "Invalid Amplicon"  
  #Find any edges that are translocations or distal by going through edge list, if none, it's non-circular
  max_path = 0
  is_complex = False
  translocation = False
  total_size = sum([s.end-s.start for s in amplicon['graph']['segments']])
  for edge in amplicon['graph']['discordant_edges']:
    if edge[0].chrom in decoys or edge[1].chrom in decoys:
        continue
    elif edge[0].chrom != edge[1].chrom:
      translocation = True
    else:
      max_path = max(abs(edge[0].start-edge[1].start),max_path)  
  if (max_path > 1e6 or translocation):
    return "Heavily-rearranged"
  if int(amplicon['amplicon']['#CoverageShiftsWithBreakpointEdges']) > 0:
    return "Non-circular"
  return "Non-circular"

def compute_av_amplification(cnvs,seg, keyword='Amplification'):
  hits = [h[0] for h in cnvs.intersection([seg])]
  hits.sort()
  pos = 0
  avAmp = 0
  for h in hits:
    start = seg.start if h.start <= seg.start else h.start
    end = seg.end if h.end >= seg.end else h.end
    pos+=end-start+1
    avAmp+=(end-start+1)*h.info[keyword] 
  avAmp+=(seg.end-seg.start+1)-pos
  avAmp/=(seg.end-seg.start+1)
  return avAmp

def read_graph_file(file, chr_prefix = True):
  input = open(file,'r')
  foo = input.next()
  segments = hg19.interval_list()
  discordant_edges = []
  concordant_edges = []
  for line in input:
    res = line.split('\t')
    if res[0] == 'sequence':
      chrom = "%s%s" % ('chr' if chr_prefix else "", res[1].split(':')[0])    
      start = min(int(res[1].split(':')[1][0:-1]),int(res[2].split(':')[1][0:-1]))
      end = max(int(res[1].split(':')[1][0:-1]),int(res[2].split(':')[1][0:-1]))    
      segments.append(hg19.interval(chrom, start, end, info={'Amplification':float(res[3])/2.}))
      segments.sort()
    if res[0] == 'discordant':
      (start, end)  = res[1].split('>')
      start = hg19.interval("%s%s-%s" % ("chr" if chr_prefix else "",start[0:-2],start.split(':')[1][0:-2]), info={'orientation':start[-1],'AvAmp':float(res[2])/2})
      end = hg19.interval("%s%s-%s" % ("chr" if chr_prefix else "",end[0:-1],end.split(':')[1][0:-1]), info={'orientation':end[-1],'AvAmp':float(res[2])/2})
      discordant_edges.append((start,end))
    if res[0] == 'concordant':
      (start, end)  = res[1].split('>')
      start = hg19.interval("%s%s-%s" % ("chr" if chr_prefix else "",start[0:-2],start.split(':')[1][0:-2]), info={'orientation':start[-1],'AvAmp':float(res[2])/2})
      end = hg19.interval("%s%s-%s" % ("chr" if chr_prefix else "",end[0:-1],end.split(':')[1][0:-1]), info={'orientation':end[-1],'AvAmp':float(res[2])/2})
      concordant_edges.append((start,end))      
  segments.sort()
  prev = len(segments)
  reduced = reduce_regions(segments)
  if len(reduced) != prev:
    print "%s is bad" % file
  input.close()
  reduced.sort()
  graph = {'segments':reduced,'discordant_edges':discordant_edges,'concordant_edges':concordant_edges}
  return graph

def reduce_regions(regions, threshold=0):
  if len(regions) == 0:
    return regions
  regions.sort()
  new_region = hg19.interval_list([copy.copy(regions.pop(0))])
  while (len(regions) > 0):
    next = copy.copy(regions.pop(0))
    if new_region[-1].intersects(next,threshold):
      new_region[-1].start = min(new_region[-1].start,next.start)
      new_region[-1].end = max(new_region[-1].end,next.end)
    else:
      new_region.append(next)  
  return new_region

def read_summary_file(summary_file, cnv_file=None):
  input = open(summary_file,'r')
  amplicons = {}
  line = next(input, None)
  if line is None:
    return None
  if int(line.split(' ')[-1].strip()) == 0:
    return amplicons
  line = next(input, None)
  if line is None:
    return None
  id = None
  for line in input:
    res = line.strip().split(' ')
    if len(res) == 2:
      continue
    if res[1] == 'AmpliconID':
      id = int(res[3])
      amplicons.setdefault(id,{})      
    else:
      amplicons[id].setdefault(res[1],res[3])  
  if len(amplicons[id].keys()) < 5:
    return None
  return amplicons

def read_cycle_file(cycle_file):
  input = open(cycle_file)
  segment_map = {}
  interval_map = {}
  cycle_map = {}
  for line in input:
    res = line.strip().split('\t')
    if res[0] == 'Interval':
      interval_map[res[1]] = hg19.interval(res[2] if res[2].find('chr') != -1 else "chr%s" % res[2],int(res[3]),int(res[4]),info={'line':res})
    if res[0] == 'Segment':
      segment_map[res[1]] = hg19.interval(res[2] if res[2].find('chr') != -1 else "chr%s" % res[2],int(res[3]),int(res[4]),info={'line':res})
    if res[0].find('Cycle') != -1:
      segments = res[0].split(';')
      cycle_name = segments[0].split('=')[-1]
      copy_count = float(segments[1].split('=')[-1])
      cycle = segments[2].split('=')[-1].split(',')
      if cycle[0][0] != '0' and len([c for c in cycle if c[0] == '0']) > 0:
        idx = cycle.index([c for c in cycle if c == '0+'][0])
        cycle = cycle[idx:] + cycle[:idx]
      cycle_map[cycle_name] = {'cycle':cycle,'copy_count':copy_count}
  input.close()
  return (segment_map, interval_map, cycle_map)   

def load_aa_result_final(out_dir, prefix, chr_prefix = True):
  summary_file = "%s/%s_summary.txt" % (out_dir, prefix)
  all_map = {}
  if os.path.exists(summary_file):    
    amps = read_summary_file(summary_file)
    if amps is None:
      return None
  else:
    return None
  for (ai, a) in amps.items():
    ints = a['Intervals'].split(',')
    cycle_file = "%s/%s_amplicon%s_cycles.txt" % (out_dir,prefix,ai)
    graph_file = "%s/%s_amplicon%s_graph.txt" % (out_dir,prefix,ai)
    (segment_map, interval_map, cycle_map) = read_cycle_file(cycle_file)
    graph = read_graph_file(graph_file, chr_prefix = chr_prefix)
    all_map.setdefault(ai,{})['segment_map'] = segment_map
    all_map.setdefault(ai,{})['interval_map'] = interval_map
    all_map.setdefault(ai,{})['cycle_map'] = cycle_map
    all_map.setdefault(ai,{})['graph'] = graph
    all_map.setdefault(ai,{})['amplicon'] = amps[ai]
    for i in interval_map:
      amp = compute_av_amplification(graph['segments'],interval_map[i])  
      interval_map[i].info['AvAmp'] = amp
      interval_map[i].info['sample'] = prefix  
      interval_map[i].info['Amplicon'] = ai      
    for s in segment_map:
      amp = compute_av_amplification(graph['segments'],segment_map[s])  
      segment_map[s].info['AvAmp'] = amp      
      segment_map[s].info['sample'] = prefix      
      segment_map[s].info['Amplicon'] = ai
    classification = classify_final_amplicon(all_map[ai], threshold_copy = 4, min_length = 10000)
    all_map.setdefault(ai,{})['classification'] = classification
    amplicon_size = sum([i.end - i.start for i in interval_map.values()])
    for i in interval_map:
      interval_map[i].info['classification'] = classification
      interval_map[i].info['amplicon_size'] = amplicon_size
    for s in segment_map:
      segment_map[s].info['classification'] = classification    
      segment_map[s].info['amplicon_size'] = amplicon_size
  all_map['amplicons'] = amps  
  return all_map

def run_classify_amplicon(amplicon_dir,prefix="AA",output_path="amplicon_classification_wodecoy.txt",ref="GRCh37"):
  print("amplicon_dir=\t" + amplicon_dir)
  print("prefix=\t" + prefix)
  print("output_path=\t" + output_path)
  print("ref=\t" + ref)
  print("global_names.REF=\t" + global_names.REF)

  now = datetime.datetime.now()
  time_stamp=now.strftime('%Y-%m-%dT%H:%M:%S') + ('-%02d' % (now.microsecond / 10000))

  if os.path.exists(output_path):
    shutil.copyfile(output_path, '%s.old%s' %(output_path,time_stamp))

  fo=open(output_path,'w')
  aa_result = load_aa_result_final(amplicon_dir, prefix)
  fo.write("amplicon_dir\tprefix\tamplicon\tclassification\n")
  if aa_result:
    for amplicon in aa_result['amplicons'].keys():
      fo.write("%s\t%s\t%d\t%s\n" % (amplicon_dir, prefix, amplicon, aa_result[amplicon]['classification']))
      print("%s\t%s\t%d\t%s\n" % (amplicon_dir, prefix, amplicon, aa_result[amplicon]['classification']))
  else:
    print("NO AA RESULT in amplicon classify")    
  fo.close()
  print "FINISHED. See " + output_path
  

if __name__ == '__main__':
  import argparse
  import shutil
  import datetime
  parser = argparse.ArgumentParser()
  parser.add_argument('-d', '--dir', type=str, default='',action="store", dest="amplicon_dir")
  parser.add_argument('-p', '--prefix', type=str, default='', action='store',dest='prefix')
  parser.add_argument('-o', '--output', type=str, default='', action='store',dest='output_path')
  parser.add_argument('-r', '--ref', type=str, default='GRCh37', action='store',dest='ref')

  args = parser.parse_args()
  print("amplicon_dir=\t" + args.amplicon_dir)
  print("prefix=\t" + args.prefix)
  print("output_path=\t" + args.output_path)
  print("ref=\t" + args.ref)
  
  print("DATA_REPO=\t" + os.environ["AA_DATA_REPO"])

  if not os.path.isdir(args.amplicon_dir):
    print args.amplicon_dir + " not exist. Check if it exists"
    sys.exit(-1)
    
  import global_names
  global_names.REF = args.ref
  
  import hg19util as hg19

  run_classify_amplicon(args.amplicon_dir,prefix=args.prefix,output_path=args.output_path,ref=args.ref)

else:
  import global_names
  global_names.REF = "GRCh37"
  import hg19util as hg19



