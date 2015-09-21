# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:45:08 2015

@author: Dave
"""
import os
import pandas as pd

def make_directories(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    os.chdir(folder_path)
    if not os.path.exists("motif_variant_list"):
        os.makedirs("motif_variant_list")
    if not os.path.exists("motif_pwm_entries"):
        os.makedirs("motif_pwm_entries")
    if not os.path.exists("logo_rawfiles"):
        os.makedirs("logo_rawfiles")
    #not needed but helps with the logic of the function
    return

def run_meme(folder,fasta_filename,thread_num=40,cluster_queue="voight_mpi",machine="interdictor",meme_run="meme_v1",multithread=True,cluster=False):
    if cluster:
        if multithread:
            job_submission_str = """bsub -a openmpi -q %s -m %s -n %d -o %s/%s.log 'meme %s -protein -oc %s/%s -nostatus -time 18000 -maxsize 60000 -mod anr -nmotifs 80 -minw 5 -maxw 255 -p %d'""" % (cluster_queue,machine,thread_num,folder,folder,fasta_filename,folder,meme_run,thread_num)
        else:
            job_submission_str = """bsub -q %s -m %s -n %d -o %s/%s.log 'meme %s -protein -oc %s/%s -nostatus -time 18000 -maxsize 60000 -mod anr -nmotifs 80 -minw 5 -maxw 255'""" % (cluster_queue,machine,thread_num,folder,folder,fasta_filename,folder,meme_run)
        os.system(job_submission_str)
    else:
        meme_str = "meme %s -protein -oc %s/%s -nostatus -time 18000 -maxsize 60000 -mod anr -nmotifs 80 -minw 5 -maxw 255" % (fasta_filename,folder,meme_run)
        os.system(meme_str)
    return

#not made for cluster yet ask paul 
def run_mast(folder,fasta_filename,pwm_filename,machine="interdictor",cluster_queue="voight_normal",mast_run="mast_v1",cluster=False):
   if cluster:
       job_submission_str = """bsub -q %s -m %s -o %s/%s.log 'mast %s %s -oc %s/%s'""" % (cluster_queue,machine,folder,folder,pwm_filename,fasta_filename,folder,mast_run)
       os.system(job_submission_str)
   else:
      mast_str = "mast %s %s -oc %s/%s" % (pwm_filename,fasta_filename,folder,mast_run)
      os.system(mast_str)
    
def write_variants(motifs,foldername):
    for motif in motifs:
        with open("%s/motif_%s_list.txt" % (foldername,motif),"w") as f:
            f.write("\n".join(list(motifs[motif])))
    return
    
def create_logo_file(kmer_count_dict,folder):
    for motif in kmer_count_dict:
        with open("%s/motif_%s_logo.txt" % (folder,motif),"w") as f:
            for row in kmer_count_dict[motif][["Seq","Count"]].iterrows():
                for i in range(int(row[1]["Count"])):
                    f.write(row[1]["Seq"]+"\n")
    return

def create_motif_table(motif_file,single_file=True):
    with open(motif_file, "r") as g:
        if single_file:
            motif_dict = {line.strip().split("\t")[0]:line.strip().split("\t")[1] for line in g}
        else:
            motif_dict = [line.strip() for line in g]
    return motif_dict

def write_intervals(intervals,output_filename):
    pd.DataFrame.from_dict(intervals,orient="index").to_csv(output_filename)
    return 
    
def write_mast_pwm(pwm_dict,motifs,freq_table,alphabet,filename,psudeocount=5):
    with open(filename,"w") as f:
        index=1
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= %s\n\n" % ("".join(alphabet)))
        f.write("Background letter frequencies\n")
        freq_table = freq_table.sort("Seq")
        for row in freq_table[["Seq","Freq"]].iterrows():
            f.write(row[1]["Seq"]+" "+"%.4f "%(row[1]["Freq"]))
            if index%9==0:
                f.write("\n")
            index+=1
        f.write("\n\n")
        for pwm in sorted(pwm_dict.keys()):
            f.write("MOTIF %s %s\n" % (pwm,motifs[pwm]))
            f.write("letter-probability matrix: alength= %d w= %d nsites= %d E= 0\n"%(len(alphabet),len(motifs[pwm]),psudeocount))
            for row in pwm_dict[pwm]:
                f.write(" ".join(map(lambda x: "%.3f" % (x),row)) + "\n")
            f.write("\n")