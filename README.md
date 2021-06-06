
# WeightedLD
The application of sequencing weighting to calculations of pairwise Linkage Disequilibrium (LD).

Given a Multiple Sequence Alignment (MSA) or Variant Call File (VCF) the program will calculate per sequence its weight using the Henikoff* methodology, which is then employed in calculating the observed and expected allele frequencies , rather than treating each sequence evenly, to produce LD metrics D, D' & r2. This can be important to account for sampling bias, population structures. We demonstrate in a preprint* that this reduces the effect of uneven sampling, as underrepresented groups of sequences will each contribute more individually than redundant,
similar sequences.

-- *our paper link goes here*
-- Henikoff, S. & Henikoff, J. G. Position-based sequence weights. _Journal of Molecular Biology_  **243**, 574–578 (1994).

## Usage
We provide two implementations of WeightedLD, one in Python and another in Rust.
Both should give exactly the same numeric answers, the Rust program however is orders of magnitude faster and will be required for analyses with greater than circa 1000 variable sites.

**Python**

The Python implementation lives in `./WeightedLD.py`. Full usage instructions can be found by running `python WeightedLD.py --help`
At its most basic usage, only the  `--input` argument is required e.g.

    python WeightedLD.py --input my_alignment.fasta 
   
**Rust**

The Rust implementation lives in `./rust/weighted_ld` and requires compilation before running, see Installation below. Full usage instruction can be found by running `weighted_ld --help`
At its most basic usage, only the `--input` and `--pair-output` arguments are required.

    rust/weighted_ld/target/release/weighted_ld --input my_alignment.fasta --pair-output my_alignment_ld.tab
     
 **!note:** the rust implementation is orders of magnitude faster than the python, and can generate very large output files with relative ease given large, complex data.
 
 **Common arguments:**

`--min-acgt` Sets a minimum fraction of A,C,G & T required for a site to be considered in LD and weighting calculations, increase to account for poor sequence coverage, default: 0.8  

`--min-variability`  The minimum (dominant) minor allele fraction for a site to be considered in LD calculations, default: 0.02

`--unweighted` Use unit weights instead of Henikoff weights

 `--r2-threshold` defines a minimum squared Pearson correlation required per pairwise calculation for it to be returned in the output, default: 0.1

## Installation

Get the WeightedLD source code:
> wget https://github.com/ojcharles/WeightedLD/archive/refs/heads/main.zip
> unzip main.zip
> cd WeightedLD-main


### Python


In order to run the program, you must have a python environment where the following packages are available:
- numpy
- biopython

A simple conda environment is provided which allows users to run the program with tested dependencies.

> conda env create -f WeightedLD_conda.yml

### Rust

To build and run the Rust implementation from source:

1. Install a Rust toolchain. We recommend using [rustup](https://rustup.rs/), though your Linux distribution's package manager may have a sufficiently up to date version for you to install instead.

2.  `cd` into `./rust/weighted_ld`

3. Run `cargo build --release`. Note that this step may take a few minutes the first time.

4. Once complete, the built executable can be found at `./target/release/weighted_ld` (`./target/release/weighted_ld.exe` on Windows.

Full usage instructions can be found by running the built executable with the `--help` argument.

### More Rust performance

There is an optional feature that can be enabled in the Rust implemetation which can *substantially* increase the calculation performance (~2-10x).

Unfortunately it relies on some features of Rust which are currently only found in Rust's nightly toolchain.

To build the Rust impl with this feature:

1. Switch to the nightly toolchain by running `rustup default nightly`.

2. Run `cargo build --release --features simd`.

3. The built executable should the be available at the same location as before.

### Even more Rust performance

The above steps all end up producing a native executable that should be portable across all x86 machines running the same OS as the build machine.

There is a way to instruct the Rust compiler to produce a binary that makes full use of the current CPU's feature set.

Doing this often yields a fair degree of extra performance (up to ~50% in our testing) at the expense of producing a binary that is not guaranteed to run on CPUs other than the one used to build it.

The steps for doing this are as follows:

1. Set the environment variable `RUSTFLAGS` to `-C target-cpu=native`

- Eg in Bash, `export RUSTFLAGS='-C target-cpu=native'`

- Eg in Powershell, `$Env:RUSTFLAGS='-C target-cpu=native'`

2. Build and run as described in one of the above sections

Remember to unset this environment variable if you would like to build a portable executable again.

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on the GitHub Issues page. For questions and other
discussions feel free to contact. [Oscar Charles -
maintainer](mailto:oscar.charles.18@ucl.ac.uk)