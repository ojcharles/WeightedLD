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
        alignment: 2d numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (a, c, g, t, -, ambiguous) is
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
    
    
def henikoff_weighting(alignment_array, var_sites, minACGT):
    # todo treat val of 4 as not ok.
    # generate arrays of nSeqs
    # for each pos
        # sum the weighting to seqs
    # then normalise

    nSeqs = alignment_array.shape[0]
    nSites = alignment_array.shape[1]

    weights = np.zeros(nSeqs)    # the output
    fracOK = np.zeros(nSeqs)     
    nSitesCounted = 0                                           # we might not count all sites in the weighting
    okBaseVals = [0,1,2,3,4]                                      # which values are ok
    
    #---------- loop over each site, and for each seq keep tally of cumulative weighting
    for iSite in var_sites:                                   # for each variable site
        #print(iSite)
        array = alignment_array[:, iSite]                      # the already converted array i.e. actg--> 01234
        unique_elements, counts_elements = np.unique(array, return_counts=True)
        okBase = np.in1d(array, okBaseVals)                     # vector of T/F for okayness ignores anythin other than actg-.
        nSitesCounted = nSitesCounted + 1
        countBase = np.zeros(5, dtype=np.int64)
        countBase[0] = np.count_nonzero(array == 0)
        countBase[1] = np.count_nonzero(array == 1)
        countBase[2] = np.count_nonzero(array == 2)
        countBase[3] = np.count_nonzero(array == 3)
        countBase[4] = np.count_nonzero(array == 4)
        tSeqs = countBase[0] + countBase[1] + countBase[2] + countBase[3] + countBase[4]
        if ((tSeqs - countBase[4]) / nSeqs) < minACGT:   #if too many missing vals then go to next
            continue

        sumSiteWeight = np.uint64(0.0) # for this site whats the total weight applied to all ok seqs?
        for iSeq in range(nSeqs):  # //todo update this so that it defulats to 0 and only checks okbases
            if okBase[iSeq]: # calculate the site contribution
                iSeq_base = array[iSeq]
                siteDenom = np.uint64(tSeqs * countBase[iSeq_base])
                siteContribution = np.float64(1.0/(siteDenom))  # key    contribution to the weight from this site, this is the henikoff
                weights[iSeq] = weights[iSeq] + siteContribution        # Key    For this seq, keep a tab of its cumulative weight over all sites
                sumSiteWeight = sumSiteWeight + siteContribution
                fracOK[iSeq] = fracOK[iSeq] + 1
            
        avgWeight =  sumSiteWeight / tSeqs # sum(weights) / n(weights) - mean(weight)
        
        # give any imbigious charcters the mean weighting score
        for iSeq in range(nSeqs):   # if not okbase, give it the average for this pos
            if not okBase[iSeq]: 
                weights[iSeq] += avgWeight
    # end of loop over sites
    
    
    #---------- normalise
    norm = weights / weights.max()
    fracOK = fracOK / nSitesCounted
    
    #---------- stdout
    print("seq\tfracOK\tWeight")
    for iSeq in range(nSeqs):
        print(str(iSeq) + "\t" + str(fracOK[iSeq]) + "\t" + str(norm[iSeq]))
        
    #end
    print("alles schon")
    return norm

    
def LD(alignment_array, var_sites, weights, minACGT):
    # generate LD metric D and R2
    print("posa\tposb\tD\tD'\tR2") # stdout headers
    iLoops = 0
    outer_loop = var_sites[0:len(var_sites)-1]
    print(var_sites)
    for i in outer_loop:
        iLoops += 1
        inner_loop = var_sites[iLoops:len(var_sites)]
        for j in inner_loop:
            #print(str(i) + "   " + str(j))
            # compare each pairwise position
            
            #----- remove any sequences with indel or illegal character
            # which index values need removing?
            i_array = alignment_array[:,i]
            remove = np.where(i_array >= 4) # which seqs contain illegal char in i
            j_array = alignment_array[:,j]
            remove = np.append(remove, np.where(j_array >= 4)) # which sequences in total contain illegals?
            remove = np.unique(remove)
            
            # remove any sequence as identified above, which contained a illegal in i or j
            i_array = np.delete(i_array, remove)
            j_array = np.delete(j_array, remove)
            tSeqs = len(i_array) # for this comparison, how many seqs are there?
            # now so long as i delete the save indexes from weights im ok
            tWeights = np.delete(weights, remove)
            
            
            
           
            
            # ----- Find and Assign Major / Minor Allele
            # first site
            unique_elements, counts_elements = np.unique(i_array, return_counts=True)
            if len(unique_elements) == 1: # is the site still variable after removing?
                continue
            major = unique_elements[counts_elements.argmax()] # which value is max
            i_array = np.where(i_array == major, 1, 0) # if is major value then 1 else 0
            
            # rather than count, here we need to sum the weights. i.e. weights[np.which == 1]
            PA = sum(tWeights[i_array == 1]) / sum(tWeights) # div by su of weights for all seqs
            Pa = sum(tWeights[i_array == 0]) / sum(tWeights)
            
            # second site
            unique_elements, counts_elements = np.unique(j_array, return_counts=True)
            if len(unique_elements) == 1: # is the site still variable after removing?
                continue
            major = unique_elements[counts_elements.argmax()] # which is max
            j_array = np.where(j_array == major, 1, 0) # if is major then
            
            PB = sum(tWeights[j_array == 1]) / sum(tWeights)
            Pb = sum(tWeights[j_array == 0]) / sum(tWeights)
            
            
            
            
            
            # ----- observed allelle frequencies
            ld_obs = np.zeros(4) # as weighting is fractional
            for k in range(tSeqs): # for each sequence, see which bin it fits in. then rather than inc by 1 . increment by weighting
                if i_array[k] == 0 and j_array[k] == 0:     #Oab
                    ld_obs[0] = ld_obs[0] + tWeights[k]
                elif i_array[k] == 1 and j_array[k] == 1:   #OAB
                    ld_obs[3] = ld_obs[3] + tWeights[k]
                elif i_array[k] == 0 and j_array[k] == 1:   #OAb
                    ld_obs[1] = ld_obs[1] + tWeights[k]
                elif i_array[k] == 1 and j_array[k] == 0:   #OaB
                    ld_obs[2] = ld_obs[2] + tWeights[k]
                else:
                    print(k)
            ld_obs = ld_obs / sum(tWeights)
            
            
            
            
            
            # ----- predicted allele freqs if 0 LD
            PAB = PA * PB
            PAb = PA * Pb
            PaB = Pa * PB
            Pab = Pa * Pb
            
            
            
            

            # ----- Caclulate D [Linkage Disequilibrium]
            # the vector is now as in the hahn molpopgen book and can be used to generate D values
            # ld_obs is pAB, p
            tD = np.zeros(4)
            #tD[0] = abs(ld_obs[3] - PAB) # MajMaj
            #tD[1] = abs(ld_obs[0] - Pab) # minmin  
            #tD[2] = abs(ld_obs[2] - PAb) #Majmin
            #tD[3] = abs(ld_obs[1] - PaB) #minMaj
            tD[0] = PAB - ld_obs[3] # MajMaj
            tD[1] = Pab - ld_obs[0] # minmin  
            tD[2] = -1 * (PAb - ld_obs[2]) #Majmin
            tD[3] = -1 * (PaB - ld_obs[1]) #minMaj
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
            print(str(i)+"\t"+str(j)+"\t"+str(round(D, 4))+"\t"+str(round(DPrime, 4))+"\t"+str(round(R2, 4)))
                


def main(args):
    logging.info("Reading FASTA data from %s", args.fasta_input)
    alignment_array = read_fasta(args.fasta_input)
    logging.info("Finished reading data. Sequence count: %s, Sequence length %s", *alignment.shape)
    
    logging.info("Computing sites of iterest (min_acgt=%s)", args.min_acgt)
    var_sites = compute_variable_sites(alignment_array, args.min_acgt)
    var_sites = np.where(var_sites)[0]
    logging.info("Found %s sites of interest", len(var_sites))
    
    logging.info("Computing Henikoff weights for each sequence")
    weightsHK = henikoff_weighting(alignment_array, var_sites, args.min_acgt)

    weights1 =  np.zeros(a.nSeqs, dtype=np.uint16)
    weights1[weights1 == 0] = 1
    weights = weights1
    
    logging.info("Computing the LD parameters")
    ld = LD(alignment_array, var_sites, weights1, args.min_acgt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WeightedLD computation tool")
    parser.add_argument("--fasta-input", type=Path, help="The source file to load", required=True)
    parser.add_argument("--min-acgt", type=float, default=0.8,
        help="Minimum fractions of ACTG at a given site for the site to be included in calculation.\
            increase this to remove more noise say 0.5")
    
    args = parser.parse_args()
    main(args)