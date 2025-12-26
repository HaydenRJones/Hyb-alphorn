# Hyb-alphorn
**A longread pipeline for homoeolog recovery, and more!**
Heavily inspired by [Hybpiper](https://github.com/mossmatters/HybPiper), but designed for use with longread (ONT / PacBio) data.

Longread based assembilies enable to recovery of whole gene sequences, including introns - unlike typical hybseq aproaches which favor the recovery of exons only. 
This tool also considers cases where there are multiple copies of genes due to allo- / autopolyploidy and tried to split them using SNPs.
By doing this we can potentially avoid chimeric gene assembilies that may result from short read only assembilies of multi-copy loci.

For a more complete discussion of the methods used in this tool, and the advantages of this approach, see (BIORXIV LINK!)[XXX]

## Dependencies
- Linux/WSL2 or macOS[^*]
[^*]: Untested!
- [whatshap](https://www.biorxiv.org/lookup/doi/10.1101/085050)
- [freebayes](https://arxiv.org/abs/1207.3907)
- [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778?login=false)
- [samtools](https://pubmed.ncbi.nlm.nih.gov/33590861/)
  
- Python >= 3.11
  - pandas
  - numpy
  - mappy
  - argparse
  - #seaborn
  - #matplotlib

## Setup
Installation and setup:
```
conda create -n alphorn
conda activate alphorn
conda install python pandas numpy mappy whatshap freebayes minimap2 samtools
git clone https://github.com/HaydenRJones/Hyb-alphorn.git 
cd Hyb-alphorn
python Hyb-alphorn.py -s SAMPLE.YAML -r REFERENCE.FASTA
```

## Usage

### Command line arguments
- `-v, --version` print version info and quit
- `-s, --sample_file` A yaml file containing the names of samples to runs, and corresponding .fastq files.
See the examples for specific format
- `-r, --reference` Target loci reference list in fasta format
- `-o, --output` Directory to the save data to.
Defaults to ./output/
- `-t, --threads` Total number of threads to use.
This is split between samples when -p is set.
Defaults to 1

(optional arguments)
- `--skip_alignment`
- `--skip_assembly`
- `--skip_phasing`

### Example sample file
Sample data is formated as a .yaml file, with the following structure:

id - "sample name"

data - "full or relative path to fastq reads for this sample"

phase_n - "number of copies to phase into"

Example file:
```
-   id: 'sample_1'
    data: './data/example/sample_1.fastq'
    phase_n: 2
-   id: 'sample_2'
    data: './data/example/sample_2.fq.gz'
    phase_n: 3
...
-   id: 'sample_n'
    data: './data/example/sample_n.fq'
    phase_n: 4

```
