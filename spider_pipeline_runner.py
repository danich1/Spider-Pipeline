import spider_utils as SPU
import spider_io_handler as SPIO
import argparse
import os
import re
from collections import defaultdict
import pdb
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

def main():
    alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    machine="interdictor"
    queue="voight_mpi"
    threads=40
    parser = argparse.ArgumentParser(prog="Spider Pipeline",description="This is the spider genome pipeline. The commands are described below.")
    parser.add_argument("-mk","--create_dirs",help="Use this to create all the necessary directories.",nargs=1,metavar="Directory_Name")
    parser.add_argument("-rme","--run_meme",help="Use this to submit meme jobs onto the cluster.",nargs=2,metavar=("Folder_Name", "AA_Seq_File"))
    parser.add_argument("-nt","--no_meme_multithread", help="Use this to run meme without it being multithreaded (Optional)", action="store_true")
    parser.add_argument("-q","--queue", help="Use this to designate which queue you want to submit the job to(Optional)",nargs=1,metavar="Queue_Name")
    parser.add_argument("-m","--machine",help="Use this to designate which machine to run (Optional)",nargs=1,metavar="Machine_Name")
    parser.add_argument("-t","--threads",help="Use this to specify how many threads you want to run (Optional)",nargs=1,metavar="Thread_Number")
    parser.add_argument("-rna","--run_name",help="Use this to specify what you want to call the output(Optional)", nargs=1,metavar="Mast/Meme_Run_Name")
    parser.add_argument("-rma","--run_mast",help="Use this to run mast.",nargs=3,metavar=("Folder_Name", "AA_Seq_File", "PWM_File"))
    parser.add_argument("-c","--cluster",help="Use this to submit jobs on the cluster.",action="store_true")
    parser.add_argument("-p","--parse_mast",help="Use this to parse MAST output.",nargs=2,metavar=("Mast_File_Path","Mast_pairwise_threshold"))
    parser.add_argument("-mc","--motif_count",help="Use this to print out the motif count.",nargs=1,metavar="Motif_Freq_Count_File_Path")
    parser.add_argument("-mco","--motif_coord",help="Use this to print out the motif coordinates.",nargs=1,metavar="Motif_Freq_Coord_File_Path")
    parser.add_argument("-vf","--variant_folder",help="Use this to print out the variant folder",nargs=1,metavar="Variant_Folder_Name")
    parser.add_argument("-pc","--paint_cassettes",help="Use this to paint cassettes from MAST output.",nargs=2,metavar=("K","Interval_Folder"))
    parser.add_argument("-pw","--pairwise", help="Use this to write out Mast's pairwise table", nargs=1,metavar="Pairwise_File_Name")
    parser.add_argument("-af","--AA_freq",help="Use this to count the Amino Acid frequencies.",nargs=3,metavar=("AA_File","Motif_File/Folder","File_Output_Path"))
    parser.add_argument("-f","--folder",help="Use this to tell the pipeline that you are using a folder with the AA_freq flag. (Optional)",nargs=1,metavar="Original_Motif_File_Path")
    parser.add_argument("-lo","--logo_output",help="Use this to specify the folder for the raw logo files to be created.",nargs=1,metavar="Logo_File")    
    args = parser.parse_args()

    if args.create_dirs:
        SPIO.make_directories(args.create_dirs[0])
        
    if args.machine:
        machine = args.machine[0]
        
    if args.threads:
        threads = int(args.threads[0])

    if args.run_meme:
        run_name="meme_v1"
        if args.run_name:
            run_name = args.run_name[0]
        if args.cluster:
            multithread=True
            if args.no_meme_multithread:
                if not args.queue:
                    parser.error("Error need to specify a queue if you don't want to multithread meme")
                queue = args.queue[0]
                multithread=False
            SPIO.run_meme(args.run_meme[0],args.run_meme[1],cluster_queue=queue,thread_num=threads,machine=machine,meme_run=run_name,multithread=multithread,cluster=True)
        else:
            SPIO.run_meme(args.run_meme[0], args.run_meme[1],meme_run=run_name)
        
    if args.run_mast:
        run_name = "mast_v1"
        if args.run_name:
            run_name = args.run_name[0]
        if args.cluster:
            SPIO.run_mast(args.run_mast[0],args.run_mast[1],args.run_mast[2],machine=machine,cluster_queue=queue,mast_run=run_name,cluster=True)
        else:
            SPIO.run_mast(args.run_mast[0],args.run_mast[1],args.run_mast[2],mast_run=run_name)

    if args.parse_mast:
        motif_map,scaffold_motif_map = SPU.grab_motifs(args.parse_mast[0],int(args.parse_mast[1]))

        if args.motif_count:
            SPIO.write_motif_freq(motif_map,scaffold_motif_map,args.motif_count[0],'counts')
        if args.motif_coord:
            SPIO.write_motif_freq(motif_map,scaffold_motif_map,args.motif_coord[0],'coord')
        if args.variant_folder:
            SPIO.write_variants(scaffold_motif_map,args.variant_folder[0])
        if args.pairwise:
            print "Currently not implemented"
            #SPIO.write_pairwise(pairwise_table,args.pairwise[0])
        if args.paint_cassettes:
            super_cassette_count = defaultdict(dict)
            cassette_count = defaultdict(dict)
            motif_group_map = defaultdict(list)
            for key in motif_map.keys():
                motif_group_map[key[0]].append(key)
            colormap = plt.get_cmap('Dark2')
            cNorm = colors.Normalize(vmin=0, vmax=len(motif_group_map.keys()))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=colormap)
            motif_color_map = {}
            epsilon = 0.01
            for key_index,key in enumerate(motif_group_map.keys()):
                motif_color_map[key] = scalarMap.to_rgba(key_index)
                scale_factor = len(motif_group_map[key]) + 1
                for motif_index, motif in enumerate(sorted(motif_group_map[key])):
                    motif_color_map[motif] = scalarMap.to_rgba(key_index,alpha=float(motif_index)/scale_factor + epsilon)
            for scaffold in scaffold_motif_map:
                print "Working on Scaffold: %s" % (scaffold)
                scaffold_motif_seq = [x[2] for x in scaffold_motif_map[scaffold]['Intervals']]
                k_mer_start_len = 2
                k_mer_end_len = 4
                (sm_conv_seqs,sm_cassette_region_tree) = SPU.find_conserved_sequences(scaffold_motif_seq,k_mer_start_len,k_mer_end_len,1,lambda x,limit: x <= limit)
                (sm_conv_seqs,sm_cassette_region_tree) = SPU.fill_gaps(scaffold_motif_seq,sm_cassette_region_tree)
                (conv_seqs,super_cassette_region_tree) = SPU.find_conserved_sequences(scaffold_motif_seq,int(args.paint_cassettes[0]),k_mer_end_len,-1,lambda x,limit: x > limit)
                (conv_seqs,super_cassette_region_tree) = SPU.find_super_cassettes(sm_cassette_region_tree,conv_seqs)
                cassette_count[scaffold] = {cassette:sm_conv_seqs[cassette].items() for cassette in sm_conv_seqs}
                if len(conv_seqs) > 0:
                    super_cassette_count[scaffold] = {sup_cassette:conv_seqs[sup_cassette].items() for sup_cassette in conv_seqs}
                    #super_cassette_colormap = {label:cassette for cassette,label in enumerate(sorted(super_cassette_map.keys(),key=len))}
                    #SPU.paint_seq(args.paint_cassettes[1],scaffold_motif_map[scaffold]["Length"],scaffold,scaffold_motif_map[scaffold]['Intervals'],motif_color_map,sm_cassette_region_tree,super_cassette_map,super_cassette_colormap)
                else:
                    if len(sm_conv_seqs) < 1:
                        print "WARNING NO CASSETTES FOUND!!"
                    print "WARNING NO SUPER CASSETTES FOUND!!!!\n"
            SPIO.write_cassette_counts(cassette_count,"%s/cassettes_count" % (args.paint_cassettes[1]))
            SPIO.write_cassette_counts(super_cassette_count,"%s/super_cassettes_count" % (args.paint_cassettes[1]))
            SPIO.write_cassette_coord(cassette_count,scaffold_motif_map,"%s/cassette_coord" %(args.paint_cassettes[1]))
            SPIO.write_cassette_coord(super_cassette_count,scaffold_motif_map,"%s/super_cassette_coord" %(args.paint_cassettes[1]))

    if args.AA_freq:
        if args.folder:
            motif_dict = {}
            motif_frequency_dict = {}
            motif_files = [filename for filename in os.listdir(args.AA_freq[1]) if not filename.startswith(".")]
            main_motif = SPIO.create_motif_table(args.folder[0])
            AA_freq_count = SPU.count_AA_frequencies(args.AA_freq[0],alphabet,0.0001)
            for motif_file in motif_files:
                motif_abbv = re.match(r"\w+_(\w+\d+)_",motif_file).groups()[0]
                motif = SPIO.create_motif_table(args.AA_freq[1]+"/"+motif_file,False)
                motif_frequency_dict[motif_abbv]=SPU.count_AA_frequencies(args.AA_freq[0],motif,0)
                motif_pwm = SPU.create_AA_pwm(motif,alphabet)
                motif_dict[motif_abbv] = motif_pwm
            SPIO.write_mast_pwm(motif_dict,main_motif,AA_freq_count,alphabet,args.AA_freq[2])
            if args.logo_output:
                SPIO.create_logo_file(motif_frequency_dict,args.logo_output[0])
        else:
            motifs = SPIO.create_motif_table(args.AA_freq[1])
            AA_freq_count = SPU.count_AA_frequencies(args.AA_freq[0],alphabet,0.0001)
            AA_pwms = {seq:SPU.create_AA_pwm([motifs[seq]],alphabet) for seq in motifs}
            SPIO.write_mast_pwm(AA_pwms,motifs,AA_freq_count,alphabet,args.AA_freq[2]) 

if __name__=="__main__":
    main()
