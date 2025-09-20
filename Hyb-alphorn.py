#          __      __  __      
#   /\ |  |__)|__|/  \|__)|\ | 
#  /--\|__|   |  |\__/| \ | \| 
                            
# 	   g		 _▒▒▒_
# 	   ' 		_\'O'/
#    c	  t	   /\__|
# 	  '	  '   /	   |
# 	    a 	 /    / \
# 	     '_ /    /   \
# 	     \_/   _/    _\
                                                                
#  v0.1 Jones & Newmarch 2025 ==============================================
#  A longread pipeline for homoeolog recovery, and more!

# TODO: file stats

import os
#import sys
#import time
import argparse
import multiprocessing

import pandas as pd
import funcs_general as funcs
import funcs_parse_asm as parse_asm

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
    prog = 'Hyb-alphorn.py',
    description = 'a pipeline for processing and phasing HYBSeq enriched nanopore longreads',
    epilog = 'v0.1 - HJ & SN 15/09/25')
    
    parser.add_argument('-v', '--version',
                    action = 'version',
                    version = 'v0.1')
    
    parser.add_argument('-s', '--sample_list',
                    help = 'A text document containing the names of samples to runs, and corresponding .fastq files.\nSee the wiki for examples.',
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
                    action = argparse.BooleanOptionalAction)

    parser.add_argument('--skip_phasing',
                    help = 'Skip the phasing step',
                    action = argparse.BooleanOptionalAction)

    
    args = parser.parse_args()
    
    samples, data = funcs.parse_samples(args.sample_list)
    loci_list = funcs.get_loci_names(args.reference)
    
    if not os.path.exists(args.output):
        #print('making output folder')
        os.makedirs(args.output)       
    
    for i in range(len(samples)):
        
        if not os.path.exists(f'{args.output}/{samples[i]}'): os.makedirs(f'{args.output}/{samples[i]}')       
        
        if not args.skip_alignment:
        
            funcs.initial_alignment(samples[i], data[i], args.reference, args.output)
            
            initial_map_stats = funcs.get_map_stats(samples[i], args.output)
            initial_map_stats = initial_map_stats.split('\n')[:-1]
            initial_map_stats = [line.split('\t') for line in initial_map_stats]
            initial_map_stats = pd.DataFrame(data = initial_map_stats, columns = ['seq_name', 'seq_len', 'mapped_reads', 'unmapped_reads'])
            
            inital_sum_mapped   = initial_map_stats['mapped_reads'].astype(int).sum()
            inital_sum_ummapped = initial_map_stats['unmapped_reads'].astype(int).sum()        
        
        if not args.skip_assembly:
            
            funcs.split_to_gene_dirs(samples[i], loci_list, args.output)
            funcs.assemble_loci(samples[i], loci_list, args.output, 8)
            
    parse_asm.parse_flye_logs(samples, loci_list, args.output, '')
    parse_asm.resolve_multi_loci(args.reference, 0.75, args.output)
    
    funcs.concat_sequences(samples, args.output)
    parse_asm.trim_assembilies(samples, 500, args.reference, args.output)
