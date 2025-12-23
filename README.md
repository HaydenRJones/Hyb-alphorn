# Hyb-alphorn
A longread pipeline for homoeolog recovery, and more!

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
```
conda create -n alphorn
conda install python pandas numpy mappy whatshap freebayes minimap2 samtools
git clone ... ; cd Hyb-alphorn
python Hyb-alphorn.py ...
```
## Usage
### Command line arguments
- `-v, --version` print version info and quit
- `-s, --sample_file` A yaml file containing the names of samples to runs, and corresponding .fastq files.
See the examples for format
- `-r, --reference` Target loci reference list in fasta format
- `-o, --output` Directory to the save data to.
Defaults to ./output/
- `-t, --threads` Total number of threads to use.\nThis is split between samples when -p is set.
Default value is 1
