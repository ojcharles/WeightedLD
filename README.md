# WeightedLD

The application of sequencing weighting to calculations of pairwise Linkage Disequilibrium (LD).

Given a Multiple Sequence Alignment (MSA) or Variant Call File (VCF) the program will calculate per sequence its weight using the Henikoff* methodology, which is then employed in calculating the observed and expected allele frequencies , rather than treating each sequence evenly,  to produce LD metrics D, D' & r2.

-- *our paper link goes here*

-- Henikoff, S. & Henikoff, J. G. Position-based sequence weights. _Journal of Molecular Biology_ **243**, 574–578 (1994).


## installation
**Get WeightedLD**
> wget https://github.com/ojcharles/WeightedLD/archive/refs/heads/main.zip
> unzip main.zip
> cd WeightedLD-main

**Python**
Install the required dependencies with conda, or manually
> conda env create -f environment.yml

**Rust**
Install rustup https://rustup.rs/
> cd rust/weighted_ld
> cargo build --release

## usage
**Python**
At its most basic usage, only the `--file` argument is required

    python WeightedLD.py --file my_alignment.fasta
    
`--min-variability` can be used to alter the number of sites for which LD is computed, default: 0.02

`--min-acgt` gives the program a way to filter out sites for both LD and weighting calculations where there is poor sequence coverage, default: 0.8

`--unweighted` calculates "vanilla" LD scores such as those in PLINK.

Bringing these together

    python WeightedLD.py --file my_alignment.fasta -- min-variability 0.1 --min-acgt 0.95 --unweighted


**Rust**

At its most basic usage, only the `--input` and `--pair-output` arguments are required.

    weighted_ld --input my_alignment.fasta --pair-output my_alignment.wld
    
Other arguments

`--min-minor` the minimum minor allele frequency
`--min-acgt` gives the program a way to filter out sites for both LD and weighting calculations where there is poor sequence coverage, default: 0.8
`--r2-threshold`  the minimum pairwise r2 value to return in the output file

**!note:** the rust implementation is orders of magnitude faster than the python, and can generate very large output files with relative ease given large, complex data.


