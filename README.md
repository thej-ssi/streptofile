# Streptofile

### Group A Streptococcus (Streptococcus pyogenes) profiling from whole genome sequencing data

Current version of Streptofile performs
- **EMM typing** based on the emm nucleotide sequence database curated by the U.S. Centers for Disease Control and Prevention (https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/)
- **Multilocus Sequence Typing** based on the *S. pyogenes scheme* curated by pubMLST (https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef)
- **Virulence gene profiling** based on 66 known virulence factors



## Setup

### Conda install

```
conda install thej-ssi::streptofile
```


### pip install
```
pip install streptofile
```
non-python dependencies that need to be installed:
 - blast


### Usage

To run emm typing, MLST and virulence gene detection on a batch of assembly files
```
streptofile -o <output_folder> *.fasta
```

To run a subset of analyses, these can be specified in a comma-separated list using the --analyses parameter
```
streptofile -o <output_folder> --analyses emm,mlst,virulence *.fasta
```
