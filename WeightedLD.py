# Oscar Charles 210404
# Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 
# using matt hahn molpopgen book
# This code calculates sequence weights using the Henikoff formula from a multiple sequuence alignment
# usage python WeightedLD.py alignment.fasta

#---------- todo
# sequence weighting - treat 4 as a not ok character
# LD - why are some values being returned infinte?

import logging
import argparse
from pathlib import Path

from Bio import AlignIO
import numpy as np
import weighted_ld_helper as helper


logging.basicConfig(
    format='[%(levelname)s] %(asctime)s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)

def read_fasta(filename: Path) -> np.ndarray:
    raw_alignment = AlignIO.read(filename, "fasta")
    return helper.multiseq_to_arr(raw_alignment)


def compute_variable_sites(alignment: np.ndarray, min_acgt: float) -> np.ndarray:
    """
    Which sites have both sufficient data and multiple non-ambiguous symbols.

    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
        min_acgt: The minimum fraction of sequences which must have A/C/G/T
            symbols at a given site

    Returns:
        A 1D boolean array of length n_sites, which is True for each site that
        meets the criteria
    """
    
    # For each site, the fraction of sequences which have a concrete symbol at that position
    concrete_fraction = (alignment < 4).sum(axis=0) / alignment.shape[0]
    
    # Whether a given site has enough data to be considered for further processing
    sufficient_data = concrete_fraction > min_acgt

    # For each site, whether any sequence has a given symbol at that position
    any_acgt = [(alignment == x).any(axis=0) for x in (0, 1, 2, 3)]
    
    # Whether each site has multiple non-ambiguous symbols
    multiple_non_ambiguous = np.sum(any_acgt, axis=0) > 1
    
    return multiple_non_ambiguous & sufficient_data
    

def henikoff_weighting(alignment: np.ndarray) -> np.ndarray:
    """
    The Henikoff weighting for each sequence in the given alignment array.

    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
            - Only sites of interest are included
    
    Returns:
        A 1D numpy array of length n_seqs, containing the Henikoff weight of each sequence.
    """
    
    n_sites = alignment.shape[1]
    n_seqs = alignment.shape[0]
    
    # Mask for non-ambiguous bases (False where the base is ambiguous)
    ok_base = (alignment != 5)
    
    # For each site, the count of sequences that have a given symbol at that site
    # Eg count_base[1, 1234] contains the count of sequences that have symbol 1 at site 1234
    count_base = np.zeros((6, n_sites))
    for base in range(6):
        count_base[base, :] = (alignment == base).sum(axis=0)
        
    # For each site, the count of sequences with a concrete symbol at that site
    t_seqs = np.sum(count_base[:5, :], axis=0)
    
    site_contribution = np.zeros(alignment.shape)
    site_contribution[ok_base] = 1 / (t_seqs * count_base[alignment, np.arange(n_sites)])[ok_base]
    
    # For each ambiguous base, fill in the average contribution across all
    # sequences with concrete bases for that site
    site_contribution[~ok_base] = 0
    site_average_weight = site_contribution.sum(axis=0) / t_seqs
    site_contribution[~ok_base] = np.broadcast_to(site_average_weight, site_contribution.shape)[~ok_base]
    
    # For each sequence, the sum of contributions from each site
    weights = site_contribution.sum(axis=1)

    # Normalize such that the largest weight is exactly 1
    return weights / weights.max()

    
def ld(alignment, weights, site_map):
    """
    Generate LD metrics; D, D', and R2
    
    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
            - Only sites of interest should be included
        weights: 1D numpy array of length n_seqs, containing the relative weight of each sequence
        site_map: Some object that can map site indices in the alignment argument to meaningful site indices
    
    Returns:
        Nothing, just prints the results to stdout in a tab-separated format
    """
    
    n_seqs = alignment.shape[0]
    n_sites = alignment.shape[1]
    
    # stdout headers
    print("posa\tposb\tD\tD'\tR2")
    for first_site in range(10):
        logging.info("    Outer loop: %s/%s", first_site, n_sites)
        for second_site in range(first_site + 1, n_sites):
            # Form an array which contains all the sequences, but only the two target sites
            target_sites = alignment[:, (first_site, second_site)]

            # Remove all sequences with a bad symbol at either target site
            good_sequences = (target_sites < 4).all(axis=1)
            target_sites = target_sites[good_sequences, :]
            target_weights = weights[good_sequences]
            target_seqs = target_sites.shape[0]
            
            # For each of the sequence, the "major allele" is the symbol that occurs most frequently

            # Whether the given sequence is equal to the major symbol at the given site
            target_sites_major = np.zeros_like(target_sites, dtype=np.bool8)

            skip_site = False
            for site in (0, 1):
                unique_elements, counts = np.unique(target_sites[:, site], return_counts=True)
                if len(unique_elements) == 1:
                    # After removing the bad sequences, one or both of the
                    # sites may no longer be variable. Stop calculations here
                    # if that is the case.
                    skip_site = True
                major_symbol = unique_elements[counts.argmax()]
                target_sites_major[target_sites[:, site] == major_symbol, site] = True
            if skip_site:
                continue
            
            total_weight = target_weights.sum()
            PA, PB = np.ma.masked_array(target_weights.reshape(-1, 1).repeat(2, axis=1), ~target_sites_major).sum(axis=0) / total_weight
            Pa, Pb = np.ma.masked_array(target_weights.reshape(-1, 1).repeat(2, axis=1), target_sites_major).sum(axis=0) / total_weight
            

            # ----- predicted allele freqs if 0 LD
            PAB = PA * PB
            PAb = PA * Pb
            PaB = Pa * PB
            Pab = Pa * Pb

            
            # ----- observed allelle frequencies
            ld_obs = np.zeros(4)
            ld_obs[0] = target_weights[~target_sites_major[:, 0] & ~target_sites_major[:, 1]].sum()
            ld_obs[3] = target_weights[target_sites_major[:, 0] & target_sites_major[:, 1]].sum()
            ld_obs[1] = target_weights[~target_sites_major[:, 0] & target_sites_major[:, 1]].sum()
            ld_obs[2] = target_weights[target_sites_major[:, 0] & ~target_sites_major[:, 1]].sum()
            ld_obs = ld_obs / total_weight

            # ----- Caclulate D [Linkage Disequilibrium]
            # the vector is now as in the hahn molpopgen book and can be used to generate D values
            # ld_obs is pAB, p
            tD = np.zeros(4)
            tD[0] = PAB - ld_obs[3] # MajMaj
            tD[1] = Pab - ld_obs[0] # minmin  
            tD[2] = -1 * (PAb - ld_obs[2]) # Majmin
            tD[3] = -1 * (PaB - ld_obs[1]) # minMaj
            D = (tD[0] + tD[1] + tD[2] + tD[3]) / 4     #they should be the same anyhow

            # normalised D = D'
            if D < 0:
                denominator = max([-ld_obs[0],-ld_obs[3]])
                if denominator == 0:
                    denominator = min([-ld_obs[0],-ld_obs[3]])
            else:
                denominator = min([ld_obs[1],ld_obs[2]])
                if denominator == 0:
                    denominator = max([ld_obs[1],ld_obs[2]])  
            DPrime = D / denominator
                
            # calculate R2
            R2 = D**2 / ( PA * Pa * PB * Pb )
            
            # cat output
            print(f"{site_map[first_site]}\t{site_map[second_site]}\t{round(D, 4)}\t{round(DPrime, 4)}\t{round(R2, 4)}")


def main(args):
    logging.info("Reading FASTA data from %s", args.fasta_input)
    alignment = read_fasta(args.fasta_input)
    logging.info("Finished reading data. Sequence count: %s, Sequence length %s", *alignment.shape)
    
    logging.info("Computing sites of iterest (min_acgt=%s)", args.min_acgt)
    var_sites = compute_variable_sites(alignment, args.min_acgt)
    logging.info("Found %s sites of interest", var_sites.sum())
    
    # Trim down the alignment array to only include the sites of interest
    alignment = alignment[:, var_sites]
    
    # Maps site indices in the trimmed down array to site indices in the original alignment
    site_map = np.where(var_sites)[0]
    
    logging.info("Computing Henikoff weights for each sequence")
    weightsHK = henikoff_weighting(alignment)

    weights1 =  np.zeros(alignment.shape[0], dtype=np.uint16)
    weights1[weights1 == 0] = 1
    weights = weights1
    
    logging.info("Computing the LD parameters")
    ld(alignment, weights1, site_map)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WeightedLD computation tool")
    parser.add_argument("--fasta-input", type=Path, help="The source file to load", required=True)
    parser.add_argument("--min-acgt", type=float, default=0.8,
        help="Minimum fractions of ACTG at a given site for the site to be included in calculation.\
            increase this to remove more noise say 0.5")
    
    args = parser.parse_args()
    main(args)