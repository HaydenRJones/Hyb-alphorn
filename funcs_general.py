# funcs_general.py
# basic functions for nanoHYB
# HJ - XX/09/25

import os
import yaml
import subprocess
#import pandas as pd
#from Bio import SeqIO

def parse_samples_yaml(sample_file):  # This function is really a bit dumb. 
                                    # It will do for now but it should be rewritten at some point
    samples = []
    data    = []
    phase_n = []
    
    with open(sample_file, 'r') as file:
        sample_data = yaml.safe_load(file)
        
    for entry in sample_data: #print(sample)
        samples.append(entry.get('id'))
        data.append(entry.get('data'))
        phase_n.append(entry.get('phase_n'))
    
    sample_dict = dict(zip(samples, data))
    phase_dict  = dict(zip(samples, phase_n))
    
    return(sample_dict, phase_dict)

def parse_samples(sample_file):
    
    file  = open(sample_file, 'r')
    lines = file.readlines()
    
    samples = []
    data    = []
    
    for line in lines:
    
        split_line = line.split(':')
        samples.append(split_line[0].rstrip()) # Took me forever to catch these damn newlines!
        data.append(split_line[1].rstrip())    # :(((((
    
    sample_dict = dict(zip(samples, data))
    
    return(sample_dict)
    #return(samples, data)
    
def get_loci_names(ref_file):
    
    file  = open(ref_file, 'r')
    lines = file.readlines()
    
    loci  = [line[1:-1] for line in lines if line.startswith('>')]
    
    return(loci)
    
def make_alignment(sample, sample_data, ref_file, out_path):
    
    #minimap_command = f'minimap2 -ax map-ont {ref_file} {sample_data} | samtools view -b -F 4 | samtools sort -o {out_dir}/{sample_name}.bam --write-index -'
    minimap_command = f'minimap2 -ax map-ont {ref_file} {sample_data} | samtools sort -o {out_path} --write-index -'
    subprocess.run(minimap_command, shell = True, stdout = subprocess.PIPE, check = True)
    
    return()
    
def get_map_stats(bam_path):
    
    idxstats_command = f'samtools idxstats {bam_path}'
    stats = subprocess.run(idxstats_command, shell = True, capture_output = True, text = True)
    
    return(stats.stdout)

def concat_sequences(samples, out_dir):
    
    for sample in samples:
        concat_command = f'cat {out_dir}/{sample}/*/resolved.fasta >> {out_dir}/{sample}/{sample}_assembly.fasta'
        subprocess.run(concat_command, shell = True)
    
    return()

# TODO
# def parallel_samples():
#     return()
