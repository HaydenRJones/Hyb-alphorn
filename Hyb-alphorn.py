#          __      __  __      
#   /\ |  |__)|__|/  \|__)|\ | 
#  /--\|__|   |  |\__/| \ | \| 
                            
# 	   g		     _▒▒▒_
# 	   ' 		_\'O'/
#    c	  t	   /\__|
# 	  '	  '   /	   |
# 	    a 	 /    / \
# 	     '_ /    /   \
# 	     \_/   _/    _\
                                                               
#  v0.1 Jones & Newmarch 2025 ==============================================
#  A longread pipeline for homoeolog recovery, and more!

# TODO: logging
    
import os
#import sys
#import time
import argparse
#import multiprocessing

import pandas as pd
import funcs_phase as phase
import funcs_general as funcs
import funcs_parse_asm as assembly

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
    prog = 'Hyb-alphorn.py',
    description = 'a pipeline for processing and phasing HYBSeq enriched nanopore longreads',
    epilog = 'v0.1 - HJ & SN 31/12/25')
    
    parser.add_argument('-v', '--version',
                    action = 'version',
                    version = 'Hyb-alphorn : v0.1 - HJ & SN 31/12/25')
    
    parser.add_argument('-s', '--sample_file',
                    help = 'A yaml file containing the names of samples to runs, and corresponding .fastq files.\nSee the wiki for examples.',
                    required = True)
    
    parser.add_argument('-r', '--reference',
                    help = 'Target loci reference list in fasta format.',
                    required = True)

    parser.add_argument('-o', '--output', 
                    help = 'Directory to the save data to.',
                    default = f'{os.getcwd()}/output')
    
    parser.add_argument('-t', '--threads', 
                    help = 'Total number of threads to use.\nThis is split between samples when -p is set\nDefault value is 1',
                    type = int,
                    default = 1)    
    
    # parser.add_argument('-p', '--parallel_samples', 
    #                 help = 'Number of samples to procsess at once time.\nDefault value is 1',
    #                 type = int,
    #                 default = 1)   
    
    parser.add_argument('--skip_alignment',
                    help = 'Skip the alignment step',
                    action = argparse.BooleanOptionalAction,
                    default = False)
    
    parser.add_argument('--skip_assembly',
                    help = 'Skip the assembly step',
                    action = argparse.BooleanOptionalAction,
                    default = False)

    parser.add_argument('--skip_phasing',
                    help = 'Skip the phasing step',
                    action = argparse.BooleanOptionalAction,
                    default = False)

    
    args = parser.parse_args()
    
    samples, phase_info = funcs.parse_samples_yaml(args.sample_file)
    loci_list = funcs.get_loci_names(args.reference)
    
    initial_sample_stats = []
    full_sample_stats = []
    
    if args.output[-1] == '/': args.output = args.output[0:-1] # remove any trailing slashes from the output
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)       
    
    for sample, data in samples.items():
        
# MAIN FUNCTION WHERE WE CAN RUN PARALELL SAMPLES?
        
        initial_sample_stats.append(funcs.run_sample_alingment(sample, data, args.output, args.reference, loci_list, args.skip_alignment, args.skip_assembly, 'initial', args.threads))
            
    assembly.parse_flye_logs(samples, loci_list, args.output, '')
    assembly.resolve_multi_loci(args.reference, 0.75, args.output)
    
    funcs.concat_sequences(samples, args.output)
    assembly.trim_assembilies(samples, 500, args.reference, args.output) # TODO : add option to select overhang

    if not args.skip_phasing:
        
        for sample, data in samples.items():
            
            phase_n = phase_info.get(f'{sample}')
            
            full_sample_stats.append(funcs.run_sample_alingment(sample, data, args.output, f'{args.output}/{sample}/{sample}_trim.fasta', loci_list, False, True, 'full', args.threads))
            
# TODO : replace the next bit with a single function
            phase.call_snps(sample, phase_n, args.threads, args.output)
            phase.phase_snps(sample, phase_n, args.threads, args.output)
            phase.haplotag(sample, phase_n, args.threads, args.output)
            
    initial_sample_stats_df = pd.DataFrame(initial_sample_stats, columns = ['sample', 
                                           'initial_mapped_reads', 
                                           'initial_unamapped_reads', 
                                           'initial_percent_mapped',
                                           'initial_count_loci'])
    full_sample_stats_df = pd.DataFrame(full_sample_stats, columns = ['sample', 
                                           'full_mapped_reads', 
                                           'full_unamapped_reads',
                                           'full_percent_mapped',
                                           'full_count_loci'])
    
    initial_sample_stats_df.to_csv(f'{args.output}/inital_alignment_stats.tsv', sep = '\t', index = False)
    full_sample_stats_df.to_csv(f'{args.output}/final_alignment_stats.tsv', sep = '\t', index = False)
