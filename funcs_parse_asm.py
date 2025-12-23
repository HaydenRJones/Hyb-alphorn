# funcs_parse_asm.py
# originally flyeLogParser.py, resolveMultiContigs.py, and trimAssemblies.py
# parse flye output and write single contig loci to a file
# Look at multi assembly loci and select one to use
# Trim assemblies down down overhang from the reference
# HJ - XX/09/25

import os
import subprocess
import numpy as np
import mappy as mp
import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as pl
#import matplotlib.colors as mcolours

from Bio import SeqIO

###############################################################################

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

def get_asm_info(loci, sample, out_dir):
    
    try:
        
        file = pd.read_csv(f'{out_dir}/{sample}/{loci}/assembly_info.txt', sep = '\t', header = 0)
        n_contigs = len(file)
    
    except:
        
        n_contigs = 0
    
    return(n_contigs)

def parse_flye_logs(sample_list, loci_list, out_dir, keep_singletons): # TODO: keep single loci across samples, plot paralog heatmap like HP?
    
    output_file = np.zeros((len(loci_list), len(sample_list)), dtype = 'int')
    output_file = pd.DataFrame(output_file, columns = sample_list, index = loci_list)
    
    for sample in sample_list:
        
        for loci in loci_list:
            
            output_file.loc[loci, sample] = get_asm_info(loci, sample, out_dir)
    
        one_contig_loci = output_file[sample]
        one_contig_loci = output_file[(one_contig_loci == 1)]
        one_contig_loci = one_contig_loci.index.values.tolist()
        
        seq_list = []
        
        for loci in one_contig_loci:
            
            sequence = list(SeqIO.parse(f'{out_dir}/{sample}/{loci}/assembly.fasta', 'fasta'))
            seq_list.append([loci, f'{sequence[0].seq}'])
            
            file = open(f'{out_dir}/{sample}/{sample}_assembly.fasta', 'w')
            for i in range(len(seq_list)):
                #file.write(f'>{sample}_{seq_list[i][0]}\n{seq_list[i][1]}\n')
                file.write(f'>{seq_list[i][0]}\n{seq_list[i][1]}\n')
    
    output_file.to_csv(f'{out_dir}/assembly_info.tsv', sep = '\t')
    output_file[~output_file.eq(1).all(axis = 1)].to_csv(f'{out_dir}/multi_assembly_loci.tsv', sep = '\t')
    
    return()

###############################################################################

def make_single_alignment(loci, sample, ref_list, out_dir):
    
    index = mp.Aligner(f'{out_dir}/{sample}/{loci}/assembly.fasta', preset = 'splice')
    
    hit_temp    = []
    cigar_array = []
    qer_length  = 0
    
    for loci_name, seq in ref_list:
        
        for hit in index.map(seq):
            
            if loci_name == loci and qer_length == 0:
                qer_length = len(seq)
            
            # Manually add in clipping cigar strings infered using the start and end of the query mapping.
            # It looks like there might be a bug in mappy as this should be reported but isn't. 
            # If the bug is fixed and clipping is reported properly !this will break!
            cg = hit.cigar
            cg.insert(0, [hit.q_st, 4])
            cg.append([qer_length - hit.q_en, 4])
            
            cigar_array.append(cg)
            
            hit_temp.append([hit.ctg,
                             loci_name,
                             hit.r_st, 
                             hit.r_en, 
                             hit.r_en - hit.r_st, 
                             hit.q_st, 
                             hit.q_en, 
                             hit.q_en - hit.q_st, 
                             hit.NM, 
                             (hit.q_en - hit.q_st) / len(seq), 
                             hit.mapq,
                             cg])
    
    hit_table = pd.DataFrame(hit_temp, columns = ['CONTIG', 
                                                  'QER_NAME', 
                                                  'REF_START', 
                                                  'REF_END' , 
                                                  'MAP_LEN', 
                                                  'QER_START', 
                                                  'QER_END', 
                                                  'QER_MAP', 
                                                  'QER_MISMATCH',  
                                                  'LEN_MAP', 
                                                  'MAPQ',
                                                  'CIGAR'])
            
    return(hit_table)

def find_contam_contigs(hit_table, contigs, loci):
    
    bad_contigs = []
    hits = list(hit_table['CONTIG'])
    unmapped_contigs = set(contigs) - set(hits)
    
    for contig in unmapped_contigs:
        bad_contigs.append([contig, 'unmapped'])
        
    for contig in hit_table['QER_NAME']:
        if contig != loci:
            bad_contigs.append([contig, 'wrong_loci'])
            
    
    return()

def find_partial_contigs(hit_table, length_thresh):
    
    partial_contigs = []
    
    for i in range(len(hit_table)):
        if hit_table.iloc[i]['LEN_MAP'] <= length_thresh:
            partial_contigs.append(hit_table.iloc[i]['CONTIG'])
    
    return()

# TODO: 
# def plot_cigar(loci, cigar_array, contig_list, query_length):
#     return()

# def plot_paralogs():
#     return()

def resolve_multi_loci(ref_file, length_thresh, out_dir):
    
    multi_contig_list = pd.read_csv(f'{out_dir}/multi_assembly_loci.tsv', sep = '\t', index_col = 0)
    sample_list       = list(multi_contig_list.columns.values)
    loci_list         = list(multi_contig_list.index) 
    
    if loci_list:
        
        ref_list = []
        for record in SeqIO.parse(ref_file, 'fasta'):
            ref_list.append([record.id, str(record.seq)])
        
        for sample in sample_list:
            
            sample_contigs = pd.DataFrame(columns = ['LOCI', 'CONTIG', 'OUTCOME'])
            
            for loci in loci_list:
                if multi_contig_list[sample][loci] and multi_contig_list[sample][loci] > 1:
                    
                    assembly_info = pd.read_csv(f'{out_dir}/{sample}/{loci}/assembly_info.txt', sep = '\t', header = 0)
                    contig_list = list(assembly_info['#seq_name'])
                        
                    contig_seqs = []
                    for record in SeqIO.parse(f'{out_dir}/{sample}/{loci}/assembly.fasta', 'fasta'):
                        contig_seqs.append([record.id, str(record.seq)])
                    contig_seq = pd.DataFrame(contig_seqs, columns = ['CONTIG', 'SEQ'])
                    
                    loci_contigs = []
                    for contig in contig_list:
                        loci_contigs.append([loci, contig, 'rejected'])
                    loci_contigs = pd.DataFrame(loci_contigs, columns = ['LOCI', 'CONTIG', 'OUTCOME'])
                    
                    hit_table = make_single_alignment(loci, sample, ref_list, out_dir)
                    for contig in list(hit_table['CONTIG']):
                        loci_contigs.loc[loci_contigs['CONTIG'] == contig, 'OUTCOME'] = 'mapped'
                    
                    contam_contigs = find_contam_contigs(hit_table, contig_list, loci) # also catches unmapped
                    for contig in contam_contigs:
                        loci_contigs.loc[loci_contigs['CONTIG'] == contig[0], 'OUTCOME'] = contig[1]
                        
                    partial_contigs = find_partial_contigs(hit_table, length_thresh)
                    for contig in partial_contigs:
                        loci_contigs.loc[loci_contigs['CONTIG'] == contig, 'OUTCOME'] = 'partial'
                        
                    remaining_contigs = loci_contigs.loc[loci_contigs['OUTCOME'] == 'mapped']    
                    if len(remaining_contigs):
                        hit_table_mapped = hit_table.loc[hit_table['CONTIG'].isin(list(remaining_contigs['CONTIG']))]
                        
                        selected_hit = hit_table_mapped.iloc[hit_table_mapped['MAPQ'].idxmax()]        
                        loci_contigs.loc[loci_contigs['CONTIG'] == selected_hit['CONTIG'], 'OUTCOME'] = 'selected'
                        keep_seq = contig_seq.loc[contig_seq['CONTIG'] == selected_hit['CONTIG']]
                        
                        file = open(f'{out_dir}/{sample}/{loci}/resolved.fasta', 'w')
                        #file.write(f">{sample}_{loci}\n{list(keep_seq['SEQ'][0])}\n")
                        file.write(f">{loci}\n{keep_seq.iloc[0]['SEQ']}\n")
                    
                    sample_contigs = pd.concat([sample_contigs, loci_contigs])
                
                sample_contigs.to_csv(f'{out_dir}/{sample}/{sample}_multi_assembly_info.tsv', sep = '\t', index = False)
                
    return()
    
###############################################################################

def trim_assembilies(samples, trim_length, reference, out_dir):
    
    for sample in samples:
        
        hit_table = pd.DataFrame(columns = ['contig_name', 
                                     'mapped_name', 
                                     'map_start', 
                                     'map_end', 
                                     'read_length',
                                     'contig_length',
                                     'contig_span', 
                                     'relative_span'])
        
        index = mp.Aligner(f'{out_dir}/{sample}/{sample}_assembly.fasta', preset = 'splice')
        
        for name, seq, qual in mp.fastx_read(reference):
            
            for hit in index.map(seq):
            
                read_length = len(seq)
                read_name   = name
                mapped_name = hit.ctg
                
                #flags = []
                
                if name == mapped_name:
                    map_s, map_e = hit.r_st, hit.r_en
                    span = map_e - map_s
                    
                    hit_table.loc[len(hit_table.index)] = [mapped_name,
                                                    read_name,
                                                    map_s,
                                                    map_e,
                                                    read_length,
                                                    hit.ctg_len,
                                                    span,
                                                    span/read_length]
                    
                else:
                    
                    hit_table.loc[len(hit_table.index)] = ['NONE',
                                                    read_name,
                                                    np.nan,
                                                    np.nan,
                                                    read_length,
                                                    np.nan,
                                                    np.nan,
                                                    span/read_length]


        new_file = []
            
        for name, seq, qual in  mp.fastx_read(f'{out_dir}/{sample}/{sample}_assembly.fasta'):
            
            hit_row = hit_table.loc[hit_table['contig_name'] == name]
            
            if len(hit_row) >= 1:
                start = max(0, int(hit_row.iloc[0]['map_start']) - trim_length)
                stop  = min(int(hit_row.iloc[0]['map_end']) + trim_length, int(hit_row.iloc[0]['contig_length']))
                
                new_file.append(f'>{name}\n{seq[start:stop]}\n')
                with open(f'{out_dir}/{sample}/{sample}_trim.fasta', 'w') as file:
                    file.writelines(new_file)
                    
    return()
