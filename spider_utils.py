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

def paint_seq(folder,seq_len,output_file,motif_seqs,motif_color_map,cassette_region_tree,super_cassette_map,super_cassette_colormap):
    motif_figure = plt.figure(figsize=(22,8))
    motif_figure_legend = plt.figure(figsize=(8,10))
    graph = motif_figure.add_subplot(111)
    graph.axis([0,seq_len+10,-5,20])
    graph.minorticks_on()
    super_colormap = plt.get_cmap('ocean')
    super_cNorm = colors.Normalize(vmin=0,vmax=len(super_cassette_map.keys())+1)
    super_scalarMap=cmx.ScalarMappable(norm=super_cNorm,cmap=super_colormap)
    cassette_length = -1
    motif_set = set({})
    for index,motif in enumerate(motif_seqs):
        cassette = cassette_region_tree.get_label(index) if cassette_region_tree else None
        if cassette:
            if cassette_length > 1:
                graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),15,color=motif_color_map[motif[2][0]]))
                cassette_length = cassette_length - 1
            else:
                #beginning of the cassette
                if cassette_length == -1:
                    beginning = motif[0]
                    cassette_length = len(cassette) - 1
                    graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),15,color=motif_color_map[motif[2][0]]))
                #end of the cassette
                else:
                    graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),15,color=motif_color_map[motif[2][0]]))
                    #graph.text(beginning,17.5,cassette)
                    #graph.plot((beginning,beginning),(16,18),'k-')
                    #graph.plot((beginning,motif[1]),(17,17),'k-')
                    #graph.plot((motif[1],motif[1]),(16,18),'k-')
                    graph.add_patch(patches.Rectangle((beginning,0),(motif[1]-beginning),15,facecolor='none',linewidth=1.5))
                    cassette_length = -1
            motif_set.add(motif[2][0])
        else:
            graph.add_patch(patches.Rectangle((motif[0],0),(motif[1]-motif[0]),10,facecolor=motif_color_map[motif[2]]))
            motif_set.add(motif[2])

    for super_cassette in super_cassette_map:
        for index in super_cassette_map[super_cassette]:
            start = motif_seqs[index[0]][0]
            end = motif_seqs[index[1]][1]
            midpoint  = (start+end)/2
            graph.add_patch(patches.Arrow(midpoint,-3,(end-midpoint),0,color=super_scalarMap.to_rgba(super_cassette_colormap[super_cassette])))
            graph.add_patch(patches.Arrow(midpoint,-3,(start-midpoint),0,color=super_scalarMap.to_rgba(super_cassette_colormap[super_cassette])))
    #motif_legend = [patches.Patch(color=scalarMap.to_rgba(len(cassette_map.keys())), label="Motif")] + [patches.Patch(color=scalarMap.to_rgba(value), label="Cassette %d (%s)" % (value,label)) for label,value in sorted(cassette_map.iteritems(),key=lambda x:x[1])] + [patches.Patch(color=super_scalarMap.to_rgba(value),label="Super Cassette %d (%s)" % (value,label)) for label,value in sorted(super_cassette_colormap.items(),key=lambda x:len(x[0]))]
    motif_legend = [patches.Patch(color=motif_color_map[motif],label='Motif %s' % (motif)) for motif in motif_set] + [patches.Patch(color=super_scalarMap.to_rgba(value),label="Super Cassette %d (%s)" % (value,label)) for label,value in sorted(super_cassette_colormap.items(),key=lambda x:len(x[0]))]
    motif_figure_legend.legend(motif_legend,[patch.get_label() for patch in motif_legend])
    motif_figure.savefig("%s/%s_paint.png" % (folder,output_file))
    motif_figure_legend.savefig("%s/%s_Legend.png" % (folder,output_file))
    plt.close()
    plt.close()

def count_AA_frequencies(fasta_file,alphabet,mode="base",episilon=0.01):
    with open(fasta_file,"r") as f:
        seq_header = f.readline().strip()
        seq = f.readline().strip()
        AA_frequency_map = pd.DataFrame({base:float(len(re.findall(base,seq)))+episilon for base in alphabet}.items(),columns=["Seq","Count"])
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
    
def grab_motifs(filename,threshold):
    #delete me
    new_mapping_table = pd.read_csv("new_group_map.csv", sep="\t")
    root = ET.parse(filename).getroot()
    motif_map = {}
    scaffold_motif_map = defaultdict(dict)
    for child in root:
        if "motifs" in child.tag:
            prev_motif = ""
            for motif in child.iter("motif"):
                #delete me
                if motif.attrib["name"] in list(new_mapping_table["OLD_Type"]):
                    motif_map[motif.attrib["id"]] = motif.attrib
                    #delete me
                    motif_map[motif.attrib["id"]]["name"] = new_mapping_table[new_mapping_table["OLD_Type"] == motif.attrib["name"]]["NEW_Type"].values[0]
                #delete me
                else:
                    print motif.attrib["name"]
            #initialize the motif_1
            """
            pairwise_matrix = defaultdict(lambda:[float(0.00)]*len(motif_map.keys()))
            pairwise_matrix["motif_1"]
            #may have to change append if number of motifs increases to infinity
            for motif in child.iter("correlation"):
                if prev_motif != motif.attrib["motif_b"]:
                    index = 0
                if float(motif.attrib["value"]) > threshold:
                    pairwise_matrix[motif.attrib["motif_b"]][index] = float(motif.attrib["value"])
                else:
                    pairwise_matrix[motif.attrib["motif_b"]][index] = float(0)
                index = index + 1
                prev_motif = motif.attrib["motif_b"]
            pairwise_table = pd.DataFrame.from_dict(pairwise_matrix)
            pairwise_columns = sorted(list(pairwise_table.columns), key=lambda x: int(x[x.index("_")+1:]))
            pairwise_table = pairwise_table[pairwise_columns]
            pairwise_table = pd.DataFrame(np.array(pairwise_table) + np.array(pairwise_table.transpose()))
            pairwise_table.columns = pairwise_columns
            pairwise_table.index = pairwise_columns
            pairwise_table["Sequence"] = map(lambda x: motif_map[x]["best_f"],pairwise_columns)
            pairwise_table = pairwise_table[["Sequence"] + pairwise_columns]
            """
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
                            motif_id = hit_information.attrib["motif"]
                            #delete me
                            if motif_id not in motif_map:
                                continue
                            pos = int(hit_information.attrib["pos"])
                            #calibrate the indicies to match python's indicies
                            calibrate_pos = pos-start
                            motif_length = int(motif_map[motif_id]["width"])
                            #convert positions into python terms
                            if motif_map[motif_id]["name"] not in motif_hit_categories:
                                motif_hit_categories[motif_map[motif_id]["name"]]={"Partial":0, "Complete":0}
                            if " " in hit_information.attrib["match"]:
                                motif_hit_categories[motif_map[motif_id]["name"]]["Partial"] += 1
                            else:
                                motif_hit_categories[motif_map[motif_id]["name"]]["Complete"] += 1
                            motif_seq_intervals.append((pos,pos+motif_length-1,motif_map[motif_id]["name"]))
                            motif_variant_map[motif_map[motif_id]["name"]].add(data[calibrate_pos:calibrate_pos+motif_length])
                scaffold_motif_map[scaffold] = {"Intervals":motif_seq_intervals,"Variants":motif_variant_map,"Length":sequence_length,"Matches":motif_hit_categories}
    motif_map = {motif_map[key]["name"]:motif_map[key]["best_f"] for key in motif_map}
    return (motif_map,scaffold_motif_map)#,pairwise_table)

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
            while(index+k) < len(sequence):
                #-1 is for the overlap
                possible_cassette = ":".join(sequence[index:index+k])
                if region_tree.fits((index,index+k-1)) and possible_cassette not in seen:
                    if possible_cassette not in chunker or chunker[possible_cassette][-1][1] < index:
                        chunker[possible_cassette].append((index,index+k-1))
                index = index + 1
            k+=inc_dec_value
        cassette_table = sorted(chunker.items(),key=lambda x:len(x[1]))
        #print cassette_table
        #print sequence
        #pdb.set_trace()
        #terminate if no other cassettes are found or if all the cassettes left occur only once
        if len(cassette_table) == 0 or all(map(lambda x: len(x[1]) == 1, cassette_table)):
            return (cassette_bag,region_tree)
        cassette_to_paint = cassette_table.pop()
        if len(cassette_to_paint[1]) > 1:
            cassette_bag[cassette_to_paint[0]] = cassette_to_paint[1]
            for index in cassette_to_paint[1]:
                region_tree.insert(index,cassette_to_paint[0])
            seen.append(cassette_to_paint[0])
        k = conserve_k
        #reset window search
        chunker=defaultdict(list)

def find_super_cassettes(sm_conv_seq_tree, conv_seq):
    super_cassettes = defaultdict(list)
    super_cassette_region_tree = KD.KD_Tree()
    intervals = [val for x in sm_conv_seq_tree.to_arr() for val in range(x[0][0],x[0][1]+1,1)]
    for possible_sc in conv_seq:
        for sp_interval in conv_seq[possible_sc]:
            if all(pd.Series(range(sp_interval[0],sp_interval[1]+1,1)).isin(pd.Series(intervals))):
                super_cassettes[possible_sc].append(sp_interval)
                super_cassette_region_tree.insert(sp_interval,possible_sc)
    return (super_cassettes,super_cassette_region_tree)

def fill_gaps(sequence,region_tree):
    intervals = region_tree.to_arr()
    gap_fix_count = defaultdict(list)
    fixed_cassette_bag = defaultdict(list)
    fixed_region_tree = KD.KD_Tree()
    filled = False
    while(not(filled)):
        #scan for gaps that need to be fixed
        for index in range(len(intervals)-1):
            #change to >1 if necessary
            if (intervals[index+1][0][0]-intervals[index][0][1]) == 2:
                gap = sequence[intervals[index][0][1]+1:intervals[index+1][0][0]]
                possible_forward_str = intervals[index][1] + gap[0]
                possible_backward_str = gap[0] + intervals[index+1][1]
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
        fixed_cassette_bag[item[1]].append(tuple(item[0]))
        fixed_region_tree.insert(tuple(item[0]),item[1])
    return (fixed_cassette_bag,fixed_region_tree)