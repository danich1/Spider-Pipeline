import xml.etree.cElementTree as ET
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.patches as patches
import matplotlib.colors as colors
import KD_Tree as KD
import re
import os
import pdb

def count_AA_frequencies(fasta_file,alphabet,episilon=0.01):
    """
    Count the amino acid frequences given a fasta file

    fast_file - the fasta file name
    alphabet -  the total possible values a string can take
    episilon - the correcting term factor
    """
    with open(fasta_file,"r") as f:
        seq_header = f.readline().strip()
        seq = f.readline().strip()
        AA_frequency_map = pd.DataFrame(
            {base:float(len(re.findall(base,seq)))+episilon for base in alphabet}.items(),
            columns=["Seq","Count"]
        )
    AA_frequency_map["Freq"] = AA_frequency_map["Count"].divide(AA_frequency_map["Count"].sum())
    return AA_frequency_map
    
def create_AA_pwm(motif,alphabet):
    """
    Create a pwm matrix from the motif given

    motif - the motif string
    alphabet - the possible values a string can take
    """
    AA_matrix = np.array([[0]*len(alphabet)]*len(motif[0]),dtype=np.float64)
    AA_map = {base:index for index,base in enumerate(alphabet)}
    for sequence in motif:
        row = 0
        for base in sequence:
            AA_matrix[row][AA_map[base]]+=1
            row = row + 1
    #convert each column into frequencies
    AA_matrix = np.transpose(
        np.divide(
            np.transpose(AA_matrix),
            np.sum(AA_matrix,axis=1)
            )
        )
    AA_matrix = np.nan_to_num(AA_matrix)
    return AA_matrix
    
def grab_motifs(filename):
    """
    grabs the motifs found within a specified xml file
    
    filename - the name of the file containing motifs.
    """
    root = ET.parse(filename).getroot()
    motif_map = {}
    scaffold_motif_map = defaultdict(dict)
    for child in root:
        if "motifs" in child.tag:
            motif_idx_map = {}
            db_index = 0
            
            for motif in child.iter("motif"):
                motif_map[motif.attrib["id"]] = motif.attrib
                motif_idx_map[str(db_index)] = motif.attrib["id"]
                db_index = db_index + 1
        
        if "sequences" in child.tag:
            for seq in child.iter("sequence"):
                motif_variant_map = defaultdict(set)
                motif_seq_intervals = []
                motif_hit_categories = defaultdict(dict)
                scaffold = seq.attrib["name"]
                sequence_length = int(seq.attrib["length"])
                
                for hits in seq.iter("seg"):
                    start = int(hits.attrib["start"])
                    
                    for hit_information in hits:
                        if "data" in hit_information.tag:
                            data = re.sub(r'\n','',hit_information.text)
                        
                        if "hit" in hit_information.tag:
                            motif_id = (
                                hit_information.attrib["motif"] if "motif" in hit_information.attrib 
                                else hit_information.attrib["idx"]
                            )
                            pos = int(hit_information.attrib["pos"])
                            motif_name = (
                                motif_map[motif_id]["name"] if not motif_id.isdigit() 
                                else motif_idx_map[motif_id]
                            )
                            
                            #calibrate the indicies to match python's indicies
                            calibrate_pos = pos-start
                            motif_length = (
                                int(motif_map[motif_id]["width"]) if not motif_id.isdigit() 
                                else int(motif_map[motif_idx_map[motif_id]]["length"])
                            )
                            
                            #convert positions into python terms
                            if motif_name not in motif_hit_categories:
                                motif_hit_categories[motif_name]={"Partial":0, "Complete":0}
                           
                            if " " in hit_information.attrib["match"]:
                                motif_hit_categories[motif_name]["Partial"] += 1
                           
                            else:
                                motif_hit_categories[motif_name]["Complete"] += 1
                           
                            motif_seq_intervals.append((pos,pos+motif_length-1,motif_name))
                            motif_variant_map[motif_name].add(data[calibrate_pos:calibrate_pos+motif_length])
                
                scaffold_motif_map[scaffold] = {
                    "Intervals":motif_seq_intervals,
                    "Variants":motif_variant_map,
                    "Length":sequence_length,
                    "Matches":motif_hit_categories
                }

    motif_map = {
        motif_map[key]["name"]:motif_map[key]["best_f"] 
        for key in motif_map} if not motif_id.isdigit() else {motif_idx_map[val]:motif_map[motif_idx_map[val]]['alt'] 
        for val in motif_idx_map
    }
    return (motif_map, scaffold_motif_map) 

def find_conserved_sequences(sequence,k,limit,inc_dec_value,condition_func):
    """
    Finds the conserved cassettes

    sequence - the amino acid sequence string
    k - the maximum or minimum size of cassettes to be found
    limit - the smallest or largest cassette to be found
    inc_dec_value - the value to change the count by
    condition_func - the function that specifies if the end has been reached.
    """

    region_tree = KD.KD_Tree()
    chunker = defaultdict(dict)
    cassette_bag = defaultdict(dict)
    seen = []
    conserve_k = k

    while(True):
        while(condition_func(k,limit)):
            index = 0

            while(index+k) < len(sequence):
                possible_group = ":".join(map(lambda x: x[0:x.index("_")],sequence[index:index+k]))
                possible_cassette = ":".join(sequence[index:index+k])
                
                if region_tree.fits((index,index+k-1)) and possible_group not in seen:
                    if possible_group not in chunker:
                        chunker[possible_group]=defaultdict(list)
                        chunker[possible_group][possible_cassette].append((index,index+k-1))
                    
                    elif possible_cassette not in chunker[possible_group]:
                        chunker[possible_group][possible_cassette].append((index,index+k-1))
                    
                    elif chunker[possible_group][possible_cassette][-1][1] < index:
                        chunker[possible_group][possible_cassette].append((index,index+k-1))
                
                index = index + 1
            k+=inc_dec_value
        
        cassette_table = sorted(
            chunker.items(),
            key=lambda x: sum(map(lambda interval: len(interval),x[1].values()))
        )

        if len(cassette_table) == 0 or all(map(lambda interval: len(interval[1].values()) == 1, cassette_table)):
            return (cassette_bag, region_tree)
        
        cassette_to_paint = cassette_table.pop()
        cassette_bag[cassette_to_paint[0]] = cassette_to_paint[1]
        
        for cassette_type in cassette_to_paint[1]:
            for index in cassette_to_paint[1][cassette_type]:
                if region_tree.fits(index):
                    region_tree.insert(index, cassette_type)
        
        seen.append(cassette_to_paint[0])
        k = conserve_k

        #reset window search
        chunker=defaultdict(dict)

def find_super_cassettes(sm_conv_seq_tree, conv_seq):
    """
    Find supercassettes after finding the cassettes

    sm_conv_seq_tree - the kd tree containing the cassette region tree
    conv_seq - list of cassettes 
    """

    super_cassettes = defaultdict(dict)
    super_cassette_region_tree = KD.KD_Tree()
    forbidden_interval = [
        val[0][1] 
        if val[0][1] - val[0][0] <= 1 else val2 
        for val in sm_conv_seq_tree.to_arr() 
        for val2 in range(val[0][0]+1,val[0][1]+1)
    ]
    
    for possible_group in conv_seq:
        for possible_sc in conv_seq[possible_group]:
            allowed = True
            
            if len(conv_seq[possible_group][possible_sc]) > 1 or len(conv_seq[possible_group]) > 1:

                for sp_interval in conv_seq[possible_group][possible_sc]:
                    if sp_interval[0] in forbidden_interval or sp_interval[1]+1 in forbidden_interval:
                        allowed = False

                    if allowed:
                        cassettes = [
                        sm_conv_seq_tree.get_label(val) 
                        for val in range(sp_interval[0],sp_interval[1]+1,1) 
                        if val not in forbidden_interval
                        ]

                        if None not in cassettes:
                            possible_group_name = "::".join(map(lambda x: ":".join(re.findall(r'([\w]+)_',x)),cassettes))
                            if possible_group_name not in super_cassettes:
                                super_cassettes[possible_group_name] = defaultdict(list)

                            if super_cassette_region_tree.fits(sp_interval):
                                super_cassettes[possible_group_name]["::".join(cassettes)].append(sp_interval)
                                super_cassette_region_tree.insert(sp_interval,"::".join(cassettes))
    
    return (super_cassettes,super_cassette_region_tree)

def fill_gaps(sequence,region_tree):
    """
    This method fills the gaps between cassettes found using meme.

    sequence- the string of amino acids
    region_tree- the KD tree that contains regions of each found cassette
    """

    intervals = region_tree.to_arr()
    gap_fix_count = defaultdict(list)
    fixed_cassette_bag = defaultdict(dict)
    fixed_region_tree = KD.KD_Tree()
    filled = False
    while(not(filled)):
        
        #scan for gaps that need to be fixed
        for index in range(len(intervals)-1):
            
            #change to >1 if necessary
            if (intervals[index+1][0][0]-intervals[index][0][1]) == 2:
                gap = sequence[intervals[index][0][1]+1:intervals[index+1][0][0]]
                possible_forward_str = intervals[index][1] + ":" + gap[0]
                possible_backward_str = gap[0] + ":" +  intervals[index+1][1]
                gap_fix_count[possible_forward_str].append([index,intervals[index][0][0],intervals[index][0][1]+1])
                gap_fix_count[possible_backward_str].append([index+1,intervals[index+1][0][0]-1,intervals[index+1][0][1]])
        
        #only consider gaps that occur more than once
        gap_fix_count = {key:gap_fix_count[key] for key in gap_fix_count if len(gap_fix_count[key]) > 1}
        
        if len(gap_fix_count) == 0:
            filled = True
        
        else:
            #pick the most freq case
            fill_val =  max(gap_fix_count.items(),key=lambda x: len(x[1]))
            for fixed_interval in fill_val[1]:
                intervals[fixed_interval[0]][0][0] = fixed_interval[1]
                intervals[fixed_interval[0]][0][1] = fixed_interval[2]
                intervals[fixed_interval[0]][1] = fill_val[0]
            gap_fix_count = defaultdict(list)
    
    for item in intervals:
        fixed_motif_group = ":".join(map(lambda x: x[0:x.index("_")],item[1].split(":")))
        if fixed_motif_group not in fixed_cassette_bag:
            fixed_cassette_bag[fixed_motif_group] = defaultdict(list)
        
        fixed_cassette_bag[fixed_motif_group][item[1]].append(tuple(item[0]))
        fixed_region_tree.insert(tuple(item[0]),item[1])
    
    return (fixed_cassette_bag,fixed_region_tree)
