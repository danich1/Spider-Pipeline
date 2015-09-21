import xml.etree.cElementTree as ET
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.cm as cmx
import KD_Tree as KD
import re
import os

def paint_seq(folder,seq_file,output_file,motif_seqs,cassette_map,cassette_region_tree,super_cassette_map,super_cassette_colormap):
    with open(seq_file,"r") as f:
        f.readline()
        seq_len = len(f.readline().strip())
    motif_figure = plt.figure(figsize=(22,8))
    motif_figure_legend = plt.figure(figsize=(8,10))
    graph = motif_figure.add_subplot(111)
    graph.axis([0,seq_len+10,-5,18])
    graph.minorticks_on()
    colormap = plt.get_cmap('Dark2')
    super_colormap = plt.get_cmap('ocean')
    cNorm = colors.Normalize(vmin=0, vmax=len(cassette_map.keys())+1)
    super_cNorm = colors.Normalize(vmin=0,vmax=len(super_cassette_map.keys())+1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=colormap)
    super_scalarMap=cmx.ScalarMappable(norm=super_cNorm,cmap=super_colormap)
    for index,motif in enumerate(motif_seqs):
        cassette = cassette_region_tree.get_label(index)
        if cassette:
            graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),15,color=scalarMap.to_rgba(cassette_map[cassette])))
        else:
            graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),10,facecolor=scalarMap.to_rgba(len(cassette_map.keys()))))
    for super_cassette in super_cassette_map:
        for index in super_cassette_map[super_cassette]:
            start = motif_seqs[index[0]][0]
            end = motif_seqs[index[1]][1]
            midpoint  = (start+end)/2
            graph.add_patch(patches.Arrow(midpoint,-3,(end-midpoint),0,color=super_scalarMap.to_rgba(super_cassette_colormap[super_cassette])))
            graph.add_patch(patches.Arrow(midpoint,-3,(start-midpoint),0,color=super_scalarMap.to_rgba(super_cassette_colormap[super_cassette])))
    motif_legend = [patches.Patch(color=scalarMap.to_rgba(len(cassette_map.keys())), label="Motif")] + [patches.Patch(color=scalarMap.to_rgba(value), label="Cassette %d (%s)" % (value,label)) for label,value in cassette_map.iteritems()] + [patches.Patch(color=super_scalarMap.to_rgba(value),label="Super Cassette %d (%s)" % (value,label)) for label,value in sorted(super_cassette_colormap.items(),key=lambda x:len(x[0]))]
    motif_figure_legend.legend(motif_legend,[patch.get_label() for patch in motif_legend])
    motif_figure.savefig("%s/%s.png" % (folder,output_file))
    motif_figure_legend.savefig("%s/%s_Legend.png" % (folder,output_file))

def count_AA_frequencies(fasta_file,alphabet,mode="base"):
    with open(fasta_file,"r") as f:
        seq_header = f.readline().strip()
        seq = f.readline().strip()
        AA_frequency_map = pd.DataFrame({base:float(len(re.findall(base,seq))) for base in alphabet}.items(),columns=["Seq","Count"])
    AA_frequency_map["Freq"] = AA_frequency_map["Count"].divide(AA_frequency_map["Count"].sum())
    return AA_frequency_map
    
def create_AA_pwm(motif,alphabet):
    AA_matrix = np.array([[0]*len(alphabet)]*len(motif[0]),dtype=np.float64)
    AA_map = {base:index for index,base in enumerate(alphabet)}
    for sequence in motif:
        row = 0
        for base in sequence:
            AA_matrix[row][AA_map[base]]+=1
            row = row + 1
    #convert each column into frequencies
    AA_matrix = np.transpose(np.divide(np.transpose(AA_matrix),np.sum(AA_matrix,axis=1)))
    AA_matrix = np.nan_to_num(AA_matrix)
    return AA_matrix
    
def grab_motifs(filename):
    root = ET.parse(filename).getroot()
    motif_map = {}
    motif_variant_map = defaultdict(set)
    motif_seq_intervals = []
    for child in root:
        if "motifs" in child.tag:
            for motif in child.iter("motif"):
                motif_map[motif.attrib["id"]] = motif.attrib
        if "sequences" in child.tag:
            for hits in child.iter("seg"):
                start = int(hits.attrib["start"])
                for hit_information in hits:
                    if "data" in hit_information.tag:
                        data = re.sub(r'\n','',hit_information.text)
                    if "hit" in hit_information.tag:
                        motif_id = hit_information.attrib["motif"]
                        #motif_seq += motif_map[motif_id]["name"][0:1] 
                        pos = int(hit_information.attrib["pos"])
                        #calibrate the indicies to match python's indicies
                        calibrate_pos = pos-start
                        motif_length = int(motif_map[motif_id]["width"])
                        #convert positions into python terms
                        motif_seq_intervals.append((pos-1,pos+motif_length-1,motif_map[motif_id]["name"][0:1]))
                        motif_variant_map[motif_map[motif_id]["name"]].add(data[calibrate_pos:calibrate_pos+motif_length])
    return (motif_variant_map,motif_seq_intervals)
    
def find_conserved_sequences_old(sequence,k,limit,inc_dec_value,condition_func):
    region_tree = KD.KD_Tree()
    chunker = defaultdict(list)
    cassette_bag = defaultdict(list)
    while(condition_func(k,limit)):
        index = 0
        #scan through sequence
        while(index+k) < len(sequence)+1:
            chunker[sequence[index:index+k]].append((index,index+k-1))
            index = index + 1
        #prune the regions first to remove overlapping
        for key in chunker:
            prev = -1
            remove_indicies = []
            for val in range(len(chunker[key])):
                update = True
                if prev < 0:
                    prev = val
                else:
                    if (lambda x,y: x[1] >= y[0])(chunker[key][prev],chunker[key][val]):
                        remove_indicies.append(chunker[key][val])
                        update = False
                    if update:
                        prev = val
            for interval in remove_indicies:
                chunker[key].remove(interval)
        #after done scanning check to see if any of the casettes have been found more than once
        for cassette,indicies in sorted(chunker.items(),key=lambda x:len(x[1]),reverse=True):
            accepted_regions = []
            #if the casette has been seen more than once
            if len(indicies) > 1:
                if len(cassette_bag.values()) > 0:
                    for index in indicies:
                        if region_tree.fits(index):
                            accepted_regions.append(index)
                    if len(accepted_regions) > 1:
                        for index in accepted_regions:
                            region_tree.insert(index,cassette)
                        cassette_bag[cassette] = accepted_regions
                else:
                    for index in indicies:
                        region_tree.insert(index,cassette)
                    cassette_bag[cassette] = indicies     
        k+=inc_dec_value
    return (cassette_bag,region_tree)

def find_conserved_sequences(sequence,k,limit,inc_dec_value,condition_func):
    region_tree = KD.KD_Tree()
    chunker = defaultdict(list)
    cassette_bag = defaultdict(list)
    percentage_covered = 0
    seen = []
    conserve_k = k
    while(True):
        while(condition_func(k,limit)):
            index = 0
            while(index+k) < len(sequence)+1:
                if region_tree.fits((index,index+k-1)) and sequence[index:index+k] not in seen:
                    chunker[sequence[index:index+k]].append((index,index+k-1))
                index = index + 1
            k+=inc_dec_value
        cassette_table = sorted(chunker.items(),key=lambda x:len(x[1]))
        #print cassette_table
        if len(cassette_table) == 0:
            print "Gene Covered(%): ",percentage_covered
            return (cassette_bag,region_tree)
        cassette_to_paint = cassette_table.pop()
        prev = -1
        remove_indicies = []
        for val in range(len(cassette_to_paint[1])):
            update = True
            if prev < 0:
                prev = val
            else:
                if (lambda x,y: x[1] >= y[0])(cassette_to_paint[1][prev],cassette_to_paint[1][val]):
                    remove_indicies.append(cassette_to_paint[1][val])
                    update = False
                if update:
                    prev = val
        for interval in remove_indicies:
            cassette_to_paint[1].remove(interval)
        accepted_regions = []
        for index in cassette_to_paint[1]:
            if region_tree.fits(index):
                accepted_regions.append(index)
        cassette_bag[cassette_to_paint[0]] = accepted_regions
        for index in accepted_regions:
            region_tree.insert(index,cassette_to_paint[0])
        percentage_covered += float(len(cassette_to_paint[0]) * len(cassette_to_paint[1]))/len(sequence)
        seen.append(cassette_to_paint[0])
        k = conserve_k
        #reset window search
        chunker=defaultdict(list)

def find_super_cassettes_old(sm_conv_seq, conv_seq):
    small_table = sm_conv_seq.keys()
    small_table_sizes = set(map(len,small_table))
    big_table = conv_seq.keys()
    print small_table
    print big_table
    super_cassettes = []
    match_indicies = []
    for db_str in big_table:
        is_supercassette = True
        multiples = filter(lambda x: len(db_str) % x == 0,small_table_sizes)
        # if super cassette len has to be a multiple of len smaller cassettes
        if multiples:
            #get all the matching indicies
            for query_str in small_table:
                for index in range(0,len(db_str),len(query_str)):
                    if db_str[index:index+len(query_str)] == query_str:
                        match_indicies.append(index)
            #sort in increasing order
            match_indicies = sorted(match_indicies)
            #check if the differences between values is in multiples
            while(len(match_indicies) > 1):
                val = match_indicies[0]
                match_indicies = match_indicies[1:]
                if (match_indicies[0] - val) not in multiples:
                    is_supercassette = False
            if len(match_indicies) > 0:
                #if the differences are consistent then we found super cassette
                if is_supercassette and len(db_str) - match_indicies[0] in multiples:
                    super_cassettes.append(db_str)
            match_indicies = []
    return super_cassettes 

def find_super_cassettes(sm_conv_seq, conv_seq):
    small_table = sm_conv_seq.keys()
    small_table_sizes = set(map(len,small_table))
    big_table = conv_seq.keys()
    super_cassettes = []
    match_indicies = []
    for db_str in big_table:
        available_index = 0
        supercassette = False
        multiples = filter(lambda x: len(db_str) % x == 0,small_table_sizes)
        # if super cassette len has to be a multiple of len smaller cassettes
        if multiples:
            while(not(supercassette) and available_index != len(multiples)):
                for val in range(available_index+1):
                    #get all the matching indicies
                    for query_str in small_table:
                        if len(query_str) == multiples[val]:
                            for index in range(0,len(db_str),len(query_str)):
                                if db_str[index:index+len(query_str)] == query_str:
                                    match_indicies.append(index)
                    #check to see if the string has been matched
                    for index in range(len(match_indicies)-1):
                        if not(sum(map(lambda x: (match_indicies[index+1] - match_indicies[index]) % x == 0, multiples[0:available_index+1]))):
                            match_indicies = []
                            break
                    if match_indicies:
                        supercassette = True
                    else:
                        available_index += 1
            if supercassette:
                super_cassettes.append(db_str)
    return super_cassettes

def fill_gaps(sequence,region_tree):
    intervals = region_tree.to_arr()
    gap_fix_count = defaultdict(int)
    fixed_cassette_bag = defaultdict(list)
    fixed_region_tree = KD.KD_Tree()
    fixed_intervals = []
    #scan for gaps that need to be fixed
    for index in range(len(intervals)-1):
        if (intervals[index+1][0][0]-intervals[index][0][1]) > 1:
            gap = sequence[intervals[index][0][1]+1:intervals[index+1][0][0]]
            gap_fix_count[intervals[index][1] + gap]+= 1
            gap_fix_count[intervals[index+1][1] + gap]+= 1
    
    for index in range(len(intervals)):
        if index + 1 != len(intervals):
            if (intervals[index+1][0][0]-intervals[index][0][1]) > 1:
                gap = sequence[intervals[index][0][1]+1:intervals[index+1][0][0]]
                if gap_fix_count[intervals[index][1] + gap] < gap_fix_count[intervals[index+1][1] + gap]:
                    fixed_intervals.append(intervals[index])
                    fixed_cassette_bag[intervals[index][1]].append(intervals[index][0])
                    fixed_region_tree.insert(intervals[index][0],intervals[index][1])
                    intervals[index+1][1] = intervals[index+1][1] + gap
                    intervals[index+1][0][0] = intervals[index+1][0][0]-1
                else:
                    intervals[index][1] = intervals[index][1] + gap
                    intervals[index][0][1] = intervals[index][0][1]+1
                    fixed_intervals.append(intervals[index])
                    fixed_cassette_bag[intervals[index][1]].append(intervals[index][0])
                    fixed_region_tree.insert(intervals[index][0],intervals[index][1])
            else:
                fixed_intervals.append(intervals[index])
                fixed_cassette_bag[intervals[index][1]].append(intervals[index][0])
                fixed_region_tree.insert(intervals[index][0],intervals[index][1])
        else:
            fixed_intervals.append(intervals[index])
            fixed_cassette_bag[intervals[index][1]].append(intervals[index][0])
            fixed_region_tree.insert(intervals[index][0],intervals[index][1])
    return (fixed_cassette_bag,fixed_region_tree)