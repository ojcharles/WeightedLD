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


# Implementations

## Python

The default implementation of the WeightedLD algorithm is written in Python, and lives in `./WeightedLD.py`. Full usage instructions can be found by running `python WeightedLD.py --help`

In order to run the program, you must have a python environment where the following packages are available:
- numpy
- biopython

## Rust

A second implementation of the WeightedLD algorithm has been written in Rust, and lives in `./rust/weighted_ld`.
This implementation should give exactly the same numeric answers as the Python implemetation, while being substantially faster.
This is at the expense of the source code being a little harder to read and a little more verbose.

To build and run the Rust implementation from source:
1. Install a Rust toolchain. We recommend using [rustup](https://rustup.rs/), though your Linux distribution's package manager may have a sufficiently up to date version for you to install instead.
2. `cd` into `./rust/weighted_ld`
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