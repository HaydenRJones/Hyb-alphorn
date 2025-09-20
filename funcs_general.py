# funcs_general.py
# basic functions for nanoHYB
# HJ - XX/09/25

import os
import subprocess
#import pandas as pd
#from Bio import SeqIO

def parse_samples(sample_file):
    
    file  = open(sample_file, 'r')
    lines = file.readlines()
    
    samples = []
    data    = []
    
    for line in lines:
    
        split_line = line.split(':')
        samples.append(split_line[0].rstrip()) # Took me forever to catch these damn newlines!
        data.append(split_line[1].rstrip())    # :(((((
    
    return(samples, data)
    
def get_loci_names(ref_file):
    
    file  = open(ref_file, 'r')
    lines = file.readlines()
    
    loci  = [line[1:-1] for line in lines if line.startswith('>')]
    
    return(loci)
    
def initial_alignment(sample, sample_data, ref_file, out_dir):
    
    #minimap_command = f'minimap2 -ax map-ont {ref_file} {sample_data} | samtools view -b -F 4 | samtools sort -o {out_dir}/{sample_name}.bam --write-index -'
    minimap_command = f'minimap2 -ax map-ont {ref_file} {sample_data} | samtools sort -o {out_dir}/{sample}/{sample}_initial.bam --write-index -'
    subprocess.run(minimap_command, shell = True, stdout = subprocess.PIPE, check = True)
    
    return()
    
def get_map_stats(sample, out_dir):
    
    idxstats_command = f'samtools idxstats {out_dir}/{sample}/{sample}_initial.bam'
    stats = subprocess.run(idxstats_command, shell = True, capture_output = True, text = True)
    
    return(stats.stdout)

def split_to_gene_dirs(sample, loci_list, out_dir):
    
    for loci in loci_list:
        
        if not os.path.exists(f'{out_dir}/{sample}/{loci}'): os.makedirs(f'{out_dir}/{sample}/{loci}')  
        
        split_command = f'samtools view -b {out_dir}/{sample}/{sample}_initial.bam {loci} | samtools fastq > {out_dir}/{sample}/{loci}/{loci}.fastq'
        subprocess.run(split_command, shell = True, stdout = subprocess.PIPE)
        
    return()

def assemble_loci(sample, loci_list, out_dir, threads): # TODO: other flye options
    
    for loci in loci_list:
        flye_command = f'flye -t {threads} --iterations 3 --nano-hq {out_dir}/{sample}/{loci}/{loci}.fastq -o {out_dir}/{sample}/{loci}'
        subprocess.run(flye_command, shell = True)
        
    return()

def concat_sequences(samples, out_dir):
    
    for sample in samples:
        concat_command = f'cat {out_dir}/{sample}/*/resolved.fasta >> {out_dir}/{sample}/{sample}_assembly.fasta'
        subprocess.run(concat_command, shell = True)
    
    return()
