# funcs_plotting.py
# HJ - 31/12/25

import os
import subprocess
import numpy as np
import mappy as mp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
import matplotlib.colors as mcolours

def plot_cigar(sample, out_dir, loci, cigar_tuple, contigs, query_length):
    
    hit_map = np.zeros((len(cigar_tuple), query_length), dtype = int)

    for i in range(len(hit_map)):
    
        start_pos = 0
    
            for string in cigar_tuple[i]:
                
# If we are dealing with a clip or match add that to our hit map. In total these should add to the length of the paftol reference so simply add as is.
# In the case of a splice (infered intron) we don't really have room for that in our map, so what we do is replace the last position in the hit map
# This means or final hit map isn't strictly true as all exons are shown as 1bp shorter, but it lets us visualse things equally.                
# CIGAR flags 1 ans 2 are indels in the reference or query. These are just skipped for the sake of simplicty
# We could probably look at adding this if it's really important but for now it's too much work
# If implemented it would be basically a copy paste of the splice check. 

                if string[1] == 0:
                    hit_map[i][start_pos:start_pos + string[0]] = 1
                    start_pos += string[0]
                if string[1] == 4:
                      hit_map[i][start_pos:start_pos +string[0]] = 0
                      start_pos += string[0]
                if string[1] == 3:
                    hit_map[i][start_pos] = 2
                    start_pos += 1                      
    
    cmap = mcolours.ListedColormap(['#914343', '#28a13f', '#6d47ad'])
                
    plt.figure(dpi=300)
    sns.heatmap(hit_map, 
                vmin = 0, 
                vmax = 2, 
                cmap = cmap, 
                cbar = False,
                yticklabels = contigs).set_title(f'{sample} : {loci}')
    
    for x in  range(1, len(hit_map)):
        plt.axhline(x, c = 'black', lw = 0.5)    
    
    plt.savefig(f'{out_dir}/{sample}/{loci}/{loci}.png')
    plt.close()
    
return()

def plot_heatmap():

return()