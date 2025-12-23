# funcs_phase.py
# phase related functions for nanoHYB
# HJ - XX/09/25

import os
import subprocess

def call_snps(sample, copy_number, threads, out_dir):
    
    # index_fasta_command = f'samtools faidx {out_dir}/{sample}/{sample}_trim.fasta'
    # subprocess.run(index_fasta_command, shell = True)
    
    # Need to fix this at some point otherwise things will be super slow!!!! /bin/sh: 1: Syntax error: "(" unexpected
    # Single threaded version will work for now however
    # freebayes_command = f'freebayes-parallel <(fasta_generate_regions.py {out_dir}/{sample}/{sample}_trim.fasta.fai 10000) {threads} -f {out_dir}/{sample}/{sample}_trim.fasta -p {copy_number} {out_dir}/{sample}/{sample}_full.bam --min-mapping-quality 50 --min-base-quality 20 --min-alternate-fraction 0.2  > {out_dir}/{sample}/{sample}_{copy_number}.vcf'
    freebayes_command = f'freebayes -f {out_dir}/{sample}/{sample}_trim.fasta -p {copy_number} {out_dir}/{sample}/{sample}_full.bam --min-mapping-quality 50 --min-base-quality 20 --min-alternate-fraction 0.2  > {out_dir}/{sample}/{sample}_{copy_number}.vcf'
    
    subprocess.run(freebayes_command, shell = True, stdout = subprocess.PIPE, check = True)
    #except subprocess.CalledProcessError as e:
    #    output = e.output
    
    return()

def phase_snps(sample, copy_number, threads, out_dir):
    
    whatshap_phase_command = f'whatshap polyphase {out_dir}/{sample}/{sample}_{copy_number}.vcf {out_dir}/{sample}/{sample}_full.bam --ploidy {copy_number} -o {out_dir}/{sample}/{sample}_{copy_number}_wh.vcf --threads {threads}  --ignore-read-groups --mapping-quality 50 --distrust-genotypes --reference {out_dir}/{sample}/{sample}_trim.fasta'
    subprocess.run(whatshap_phase_command, shell = True, check = True)
    
    processes_vcf_command = f'bgzip {out_dir}/{sample}/{sample}_{copy_number}_wh.vcf; tabix -p vcf {out_dir}/{sample}/{sample}_{copy_number}_wh.vcf.gz'
    subprocess.run(processes_vcf_command, shell = True, check = True)
    
    return()

def haplotag(sample, copy_number, threads, out_dir):
    
    whatshap_tag_command = f'whatshap haplotag -o {out_dir}/{sample}/{sample}_tag.bam {out_dir}/{sample}/{sample}_{copy_number}_wh.vcf.gz {out_dir}/{sample}/{sample}_full.bam --ignore-read-groups --ploidy {copy_number} --reference {out_dir}/{sample}/{sample}_trim.fasta --output-haplotag-list {out_dir}/{sample}/{sample}_wh.tsv'
    subprocess.run(whatshap_tag_command, shell = True, check = True)

    split_list = [f'{out_dir}/{sample}/{sample}_H{n}.bam' for n in range(1, copy_number + 1)] 

    whatshap_split_command = f"whatshap split -o {' -o '.join(split_list)} {out_dir}/{sample}/{sample}_tag.bam {out_dir}/{sample}/{sample}_wh.tsv"
    subprocess.run(whatshap_split_command, shell = True, check = True)
    
    for split in split_list:
        samtools_index_command = f'samtools index {split}'
        subprocess.run(samtools_index_command, shell = True, check = True)
        
    return()

def phased_consensus(sample, copy_number, loci_list, out_dir):
    
    for loci in loci_list:
        for n in range(1, copy_number + 1):
            print()
    
    return()