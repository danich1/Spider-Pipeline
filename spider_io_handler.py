# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:45:08 2015

@author: Dave
"""
import os
import pandas as pd
from collections import Counter
import pdb

def make_directories(folder_path):
    """
    This preprocessing function creates folders the entire pipeline will be using.

    folder_path - the path to create the folders
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    os.chdir(folder_path)
    
    if not os.path.exists("motif_variant_list"):
        os.makedirs("motif_variant_list")
    
    if not os.path.exists("motif_pwm_entries"):
        os.makedirs("motif_pwm_entries")
    
    if not os.path.exists("logo_rawfiles"):
        os.makedirs("logo_rawfiles")
    return

def run_meme(
    folder,fasta_filename,
    thread_num=40,cluster_queue="voight_mpi",
    machine="interdictor",meme_run="meme_v1",
    multithread=True,cluster=False
    ):
    """
    This function runs the meme motif search algorithm.

    folder - the folder to place the results
    fasta_filename - the fast file containing the sequences
    thread_num - number of threads if multithreaded
    cluster_queue - name of the cluser to submit jobs
    machine - the host machine 
    meme_run - the version of meme
    multithreaded - flag to indicate multithreading.
    cluster - flag to indicate using the cluster.
    """
    if cluster:
        if multithread:
            job_submission_str = f"bsub -a openmpi \
            -q {cluster_queue} -m {machine} -n {thread_num} -o {folder}/{folder}.log \
            'meme {fasta_filename} -protein -oc {folder}/{meme_run} \
            -nostatus -time 18000 -maxsize 60000 \
            -mod anr -nmotifs 80 -minw 5 -maxw 255 -p {thread_num}'" 
        else:
            job_submission_str = f"bsub -q {cluster_queue} -m {machine} -n {thread_num} -o {folder}/{folder}.log\
            'meme {fasta_filename} -protein -oc {folder}/{meme_run} \
            -nostatus -time 18000 -maxsize 60000 \
            -mod anr -nmotifs 80 -minw 5 -maxw 255'"
        os.system(job_submission_str)
    else:
        meme_str = f"meme {fasta_filename} -protein -oc {folder}/{meme_run} \
        -nostatus -time 18000 -maxsize 82000 \
        -mod anr -nmotifs 80 -minw 5 -maxw 255"
        os.system(meme_str)
    return

def run_mast(
    folder,fasta_filename,pwm_filename,
    machine="interdictor",cluster_queue="voight_normal",
    mast_run="mast_v1",cluster=False
    ):
    """
    This function runs the mast motif paint program.

    folder - the folder to place the results
    fasta_filename - the fast file containing the sequences
    pwm_filename - the pwm that mast will be using
    cluster_queue - name of the cluser to submit jobs
    machine - the host machine 
    mast_run - the version of mast
    multithreaded - flag to indicate multithreading.
    cluster - flag to indicate using the cluster.
    """

    if cluster:
       job_submission_str = f"bsub -q {cluster_queue} -m {machine} -o {folder}/{folder}.log \
       'mast {pwm_filename} {fasta_filename} -oc {folder}/{mast_run}'"
       os.system(job_submission_str)
   
    else:
      mast_str = f"mast {pwm_filename} {fasta_filename} -oc {folder}/{mast_run}"
      os.system(mast_str)
    
def write_variants(scaffold_motif_map,foldername):
    """
    This function writes variant motifs into a folder

    scaffold_motif_map - dictionary that contains scaffold motifs
    foldername- the name of the folder to hold the output
    """
    if not(os.path.exists(foldername)):
        os.mkdir(foldername)
    
    for scaffold in scaffold_motif_map:
        for motif in scaffold_motif_map[scaffold]['Variants']:
            with open(f"{foldername}/{scaffold+'_variants'}/motif_{motif}_list.txt", "w") as f:
                f.write("\n".join(list(scaffold_motif_map[scaffold]['Variants'][motif])))
    return

def write_motif_freq(motif_map, scaffold_motif_map,filename,flag):
    """
    This function writes information about motifs found

    motif_map - dictionary containing possible motifs
    scaffold_motif_map - dictionary containing a scaffold to motif map
    filename - name of the output file
    flag - a variable to indicate whether or not to print the count column
    """
    with open(f"{filename}.csv", "w") as f:
        
        if flag == 'counts':
            f.write("Motif Code,Motif Seq,Scaffold,Count,Partial Count,Complete Count\n")
        
        else:
            f.write("Motif Code,Scaffold,Start,Stop\n")
        
        for scaffold in scaffold_motif_map:
            if flag == 'counts':
                motif_count = Counter([element[2] for element in scaffold_motif_map[scaffold]['Intervals']])
                matches = scaffold_motif_map[scaffold]["Matches"]
                
                for motif in sorted(motif_count.keys()):
                    f.write(
                        f"{motif},{motif_map[motif]},{scaffold},\
                        {motif_count[motif]},{matches[motif]['Partial']},\
                        {matches[motif]['Complete']}\n"
                    )
            else:
                for element in scaffold_motif_map[scaffold]['Intervals']:
                    f.write(f"{element[2]},{scaffold},{element[0]},{element[1]}\n")
    return

def create_logo_file(kmer_count_dict,folder):
    """
    This function creates sequence logos to look at.

    kmer_count_dict - a dictionary that contains kmer counts
    folder- the folder to write the logos to
    """
    for motif in kmer_count_dict:
        with open(f"{folder}/motif_{motif}_logo.txt", "w") as f:
            for row in kmer_count_dict[motif][["Seq","Count"]].iterrows():
                for i in range(int(row[1]["Count"])):
                    f.write(f"{row[1]['Seq']}\n")
    return

def create_motif_table(motif_file,single_file=True):
    """
    This function creates a table of motifs and their respective group code

    motif_file - a file that contains the motif and their respective group code.
    single_file - a flag to indicate whether or not to use a single file or using a folder
    """
    with open(motif_file, "r") as g:
        next(g)
        if single_file:
            motif_dict = {line.strip().split("\t")[0]:line.strip().split("\t")[1] for line in g}
        else:
            motif_dict = [line.strip() for line in g]
    return motif_dict

def write_cassette_counts(intervals,output_filename):
    """
    This function writes cound cassettes to a file.

    intervals - a list of cassettes and their corresponding intervals
    output_filename - name of the file to output the cassettes
    """
    with open(f"{output_filename}.csv", "w") as f:
        f.write("Cassette Group,Cassette Seq,Scaffold,Count\n")
        for scaffold in intervals:
            for cassette in intervals[scaffold]:
                for cassette_type in intervals[scaffold][cassette]:
                    f.write(f"{cassette},{cassette_type[0]},{scaffold},{len(cassette_type[1])}\n")
    return 

def write_pairwise(pairwise_table, filename):
    raise NotImplementedError
    pairwise_table.to_csv(f"{filename}.csv")
    return

def write_mast_pwm(pwm_dict,motifs,freq_table,alphabet,filename,nsites=5):
    """
    This function writes out the mast pwm data to a file

    pwm_dict - a dictionary containing the pwm and the attached motif
    motifs - a list of motifs
    freq_table - the background frequency table
    alphabet - the list of possible amino acid values
    filename - the name to write out the file
    nsites - the value for mast output
    """
    with open(filename, "w") as f:
        index=1
        f.write("MEME version 4\n\n")
        f.write(f"ALPHABET= {''.join(alphabet)}\n\n")
        f.write("Background letter frequencies\n")
        freq_table = freq_table.sort_values("Seq")
        for row in freq_table[["Seq","Freq"]].iterrows():
            f.write(row[1]["Seq"]+" "+"%f "%(row[1]["Freq"]))
            if index%9==0:
                f.write("\n")
            index+=1
        f.write("\n\n")
        for pwm in sorted(pwm_dict.keys()):
            f.write(f"MOTIF {pwm} {motifs[pwm]}\n")
            f.write(f"letter-probability matrix: alength={len(alphabet)} w={len(motifs[pwm])} nsites={nsites} E=0\n")
            for row in pwm_dict[pwm]:
                f.write(" ".join(map(lambda x: f"{x:.3f}", row)) + "\n")
            f.write("\n")

def write_cassette_coord(cassettes,motif_coord,filename):
    """
    This function writes out cassette coordinates for each cassette.

    cassettes - the total cassettes to output
    motif_coord - the dictionary containing motif coordinates
    filename - the name of the file to output each coordinate.
    """
    with open(filename+".csv","w") as f:
        f.write("Cassette Group,Cassette Seq,Scaffold,Start,Stop\n")
        for scaffold in cassettes:
            for cassette_interval in cassettes[scaffold]:
                for coord in cassettes[scaffold][cassette_interval]:
                    for interval in coord[1]:
                        start = motif_coord[scaffold]["Intervals"][interval[0]][0]
                        stop = motif_coord[scaffold]["Intervals"][interval[1]][1]
                        f.write(f"{cassette_interval},{coord[0]},{scaffold},{start},{stop}\n")
