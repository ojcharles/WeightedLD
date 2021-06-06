requirements



1 - read alignment
    type fasta [*.fasta, *.fa]


2 - convert alignment to integer matrix
    "A","C","T","G","-","OTHER"
     0   1   2   3   4     5


identify variable sites for henikoff
    okbase = count those with [0,1,2,3] 
        we want to compute "-" in henikoff, but ignore sites with too many "-", or ambiguous bases, they may infer poor sequencing
    sufficient_data =   frac_okbase > input.min_acgt_fraction  //i.e. sum(0-3) / nSeq = 0.9 and  0.9 > 0.8 so TRUE
    has_variation   =   has any level of variation for [0,1,2,3,4] i.e. 99 A and 1 -
    no cap on max variability
    return boolean of sites which have: sufficient_data & has_variation i.e. T,F,T,T,T,F,T


identify variable sites for LD  //as in henikoff, and in addition:
    // as alignments can have three+ almost equal bases at a position
    identify major_base at site
    identify most dominant minor base. i.e. [AAACCT] Major -> A, domMinor -> C
    minor_fraction = dominant_minor_count / (dominant_minor_count + major_count)
        // This should be fine, as ignoring second anf third most important base just increases minor_fraction
    has_min_variability minor_fraction >= min_variability
    return boolean of sites: sufficient_data & has_variation * has_min_variability i.e. F,F,T,T,T,F,T
    // no need for major frac check now, as we will exclude in LD. ie limit to only two types of base and throw out everyhing thats not Maj or majMin


Calculate Henikoff sequence weights
    input is the matrix alignment including 4 and 5's.
    calculate the score per sequence treating 0,1,2,3,4 as unique ok values
    keep track of mean score per site
    sequences at a site wih 5 are given the mean value for that site
    each sequence accumulates a total weight, sum of weights across all sites
    final step is normalising, such that the most unique sequence has weight of 1


calculate weighted LD
    by sequence number, not weight contribution
        Identify Major base //should be limited to 0-3, only in an extreme case of 5 equal bases with  ambigious only in acgt will this error
        Identify dominant Minor base 0-4
    Identify ambigous bases
    remove sequences from calculation that are ambiguous, or (not in Major & dominant Minor).  union of pairwise sites
    calculate predicted PA, Pa, QB, Qb from real allele fractions
    sum up observed AB, Ab, aB, ab
    calculate D, D' and R2 per site.
    if r2 above 0.05 return else skip.


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