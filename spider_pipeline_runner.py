import spider_utils as SPU
import spider_io_handler as SPIO
import argparse
import os
import re

def main():
    #if you want to make it DNA ask me
    alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    machine="interdictor"
    queue="voight_mpi"
    threads=40
    parser = argparse.ArgumentParser(prog="Spider Pipeline",description="This is the spider genome pipeline. The commands are described below.")
    parser.add_argument("-mk","--create_dirs",help="Use this to create all the necessary directories.",nargs=1,metavar="Directory_Name")
    #if run meme folder,fasta_filename,thread_num,cluster_queue="voight_mpi",machine="interdictor",meme_run="meme_v1",multithread=True
    parser.add_argument("-rme","--run_meme",help="Use this to submit meme jobs onto the cluster.",nargs=2,metavar=("Folder_Name", "AA_Seq_File"))
    parser.add_argument("-nt","--no_meme_multithread", help="Use this to run meme without it being multithreaded (Optional)", action="store_true")
    parser.add_argument("-q","--queue", help="Use this to designate which queue you want to submit the job to(Optional)",nargs=1,metavar="Queue_Name")
    parser.add_argument("-m","--machine",help="Use this to designate which machine to run (Optional)",nargs=1,metavar="Machine_Name")
    parser.add_argument("-t","--threads",help="Use this to specify how many threads you want to run (Optional)",nargs=1,metavar="Thread_Number")
    parser.add_argument("-rna","--run_name",help="Use this to specify what you want to call the output(Optional)", nargs=1,metavar="Mast/Meme_Run_Name")
    #if run mast folder,fasta_filename,pwm_filename,mast_run="mast_v1",cluster=False
    parser.add_argument("-rma","--run_mast",help="Use this to run mast.",nargs=3,metavar=("Folder_Name", "AA_Seq_File", "PWM_File"))
    parser.add_argument("-c","--cluster",help="Use this to submit jobs on the cluster.",action="store_true")
    #if this will parse the mast output
    parser.add_argument("-p","--parse_mast",help="Use this to parse MAST output.",nargs=2,metavar=("Mast_File_Path","Folder_Name"))
    parser.add_argument("-pc","--parse_mast_chunks",help="Use this to parse MAST output and calculate conserved sequences.",nargs=5,metavar=("K","Interval_Folder","Interval_Output_File","Seq_File","Figure_File"))
    parser.add_argument("-af","--AA_freq",help="Use this to count the Amino Acid frequencies.",nargs=3,metavar=("AA_File","Motif_File/Folder","File_Output_Path"))
    parser.add_argument("-f","--folder",help="Use this to tell the pipeline that you are using a folder with the AA_freq flag. (Optional)",nargs=1,metavar=("Original_Motif_File_Path"))
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
        motifs,motif_seq = SPU.grab_motifs(args.parse_mast[0])
        if len(args.parse_mast) > 1:
            SPIO.write_variants(motifs,args.parse_mast[1])
        if args.parse_mast_chunks:
            k_mer_start_len = 2
            k_mer_end_len = 7
            (sm_conv_seqs,sm_cassette_region_tree) = SPU.find_conserved_sequences("".join(map(lambda x: x[2],motif_seq)),k_mer_start_len,k_mer_end_len,1,lambda x,limit: x <= limit)
            (sm_conv_seqs,sm_cassette_region_tree) = SPU.fill_gaps("".join(map(lambda x: x[2],motif_seq)),sm_cassette_region_tree)
            (conv_seqs,casette_region_tree) = SPU.find_conserved_sequences("".join(map(lambda x: x[2],motif_seq)),int(args.parse_mast_chunks[0]),k_mer_end_len,-1,lambda x,limit: x > limit)
            casette_map = {label:casette for casette,label in enumerate(sorted(sm_conv_seqs.keys(),key=len))}
            super_cassette_map = {cassette:conv_seqs[cassette] for cassette in SPU.find_super_cassettes(sm_conv_seqs,conv_seqs)}
            SPIO.write_intervals(sm_conv_seqs,"%s/cassettes_%s" % (args.parse_mast_chunks[1],args.parse_mast_chunks[2]))
            SPIO.write_intervals(conv_seqs,"%s/possible_cassettes_%s" % (args.parse_mast_chunks[1],args.parse_mast_chunks[2]))
            SPIO.write_intervals(super_cassette_map,"%s/super_cassettes_%s" % (args.parse_mast_chunks[1],args.parse_mast_chunks[2]))
            super_cassette_colormap = {label:cassette for cassette,label in enumerate(sorted(super_cassette_map.keys(),key=len))}
            SPU.paint_seq(args.parse_mast_chunks[1],args.parse_mast_chunks[3],args.parse_mast_chunks[4],motif_seq,casette_map,sm_cassette_region_tree,super_cassette_map,super_cassette_colormap)

    #add one more arugment to get the main motif sequence
    if args.AA_freq:
        if args.folder:
            motif_dict = {}
            motif_frequency_dict = {}
            motif_files = [filename for filename in os.listdir(args.AA_freq[1]) if not filename.startswith(".")]
            main_motif = SPIO.create_motif_table(args.folder[0])
            AA_freq_count = SPU.count_AA_frequencies(args.AA_freq[0],alphabet)
            for motif_file in motif_files:
                motif_abbv = re.match(r"\w+_(\w+\d+)_",motif_file).groups()[0]
                motif = SPIO.create_motif_table(args.AA_freq[1]+"/"+motif_file,False)
                motif_frequency_dict[motif_abbv]=SPU.count_AA_frequencies(args.AA_freq[0],motif)
                motif_pwm = SPU.create_AA_pwm(motif,alphabet)
                motif_dict[motif_abbv] = motif_pwm
            SPIO.write_mast_pwm(motif_dict,main_motif,AA_freq_count,alphabet,args.AA_freq[2])
            if args.logo_output:
                SPIO.create_logo_file(motif_frequency_dict,args.logo_output[0])
        else:
            motifs = SPIO.create_motif_table(args.AA_freq[1])
            AA_freq_count = SPU.count_AA_frequencies(args.AA_freq[0],alphabet)
            AA_pwms = {seq:SPU.create_AA_pwm([motifs[seq]],alphabet) for seq in motifs}
            SPIO.write_mast_pwm(AA_pwms,motifs,AA_freq_count,alphabet,args.AA_freq[2]) 

#makes it so other files won't execute unnecessary commands
if __name__=="__main__":
    main()
