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
    Identify Major base //should be limited to 0-3, only in an extreme case of 5 equal bases with  ambigious only in acgt will this error
    Identify dominant Minor base 0-4
    Identify ambigous bases
    remove sequences from calculation that are ambiguous, or (not in Major & dominant Minor).  union of pairwise sites
    calculate predicted PA, Pa, QB, Qb from real allele fractions
    sum up observed AB, Ab, aB, ab
    calculate D, D' and R2 per site.
    if r2 above 0.05 return else skip.


