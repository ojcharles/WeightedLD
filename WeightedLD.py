# Oscar Charles 210404
# #Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 
# This code calculates sequence weights using the Henikoff formula from a multiple sequuence alignment
# usage python WeightedLD.py alignment.fasta

import sys
from Bio import AlignIO
import os
import numpy as np

os.chdir('C:\\Oscar\\OneDrive\\UCL\\code\\WeightedLD')


### modifications
minACGT = 0.5   # Minimum fractions of ACTG at a given site for the site to be included in calculation. increase this to remove more noise say 0.5
### end

# handle args
#msa = sys.argv[0]
# alignmentFile = "all_raw_best_msa_man4.fasta"
alignmentFile = "test.fasta"



class MSA:
    ### this handles the alignment
    # read alignment
    # generate allelle matrix
    # todo
    # make the numeric array small that int23, memory efficient for larger alignments
    # henikoff weighting
    # make var sites integer
    def __init__(self, alignmentFile):
        self.alignment = AlignIO.read(alignmentFile, "fasta")
        self.nSeqs = len(self.alignment)
        self.nSites = self.alignment.get_alignment_length()
        self.alignment_array = self.alignment2alignment_array(self.alignment)
        self.var_sites = self.which_pos_var(self.alignment_array, self.nSites)
        self.weighting = self.henikoff_weighting(self)
    
    def alignment2alignment_array(self, alignment):
        alignment_array = np.array([list(rec) for rec in alignment], np.dtype("U4")) # alignment as matrix
        alignment_array = alignment_array.astype('U13') #convert from bytecode to character       
        BaseToInt = {
           "A":0,
           "a":0,
           "C":1,
           "c":1,
           "G":2,
           "g":2,
           "T":3,
           "t":3,
           "-":4,
           ",":4}
        # translate character matrix to integer using dict
        alignment_array = np.vectorize(BaseToInt.get)(alignment_array)
        return alignment_array
            
    
    def which_pos_var(self, alignment_array, nSites):
        # in
            # alignment array
        # out
            # vector of variable position ignore character 4
        which_var = np.zeros(shape=(nSites,1),dtype=(np.int32))
        for pos in range(0,nSites):
            #print(pos)
            a = np.unique(alignment_array[:,pos])
            a = a[a!= 4] # remove the indels
            a = len(a)
            if a > 1:
                which_var[pos] = pos
        out = which_var[which_var != 0] 
        return out

    
    def henikoff_weighting(self, nSeqs):
        # todo
        return np.zeros(shape=(self.nSeqs,1))
        
    
    
    def LD_D(self, alignment_array, var_sites):
        # generate LD metric D
        return 1
        
a = MSA(alignmentFile)







outer_loop = a.var_sites[0:len(a.var_sites)-1]
inner_loop = a.var_sites[1:len(a.var_sites)]
for i in outer_loop:
    for j in inner_loop:
        print(str(i) + "   " + str(j))
        # compare each pairwise position
        
        # first site
        # find major allele
        i_array = a.alignment_array[:,i]
        unique_elements, counts_elements = np.unique(i_array, return_counts=True)
        major = unique_elements[counts_elements.argmax()] # which is max
        i_array = np.where(i_array== major, 1, 0) # if is major then
                
        
        # second site
        # find major allelle
        j_array = a.alignment_array[:,j]
        unique_elements, counts_elements = np.unique(i_array, return_counts=True)
        major = unique_elements[counts_elements.argmax()] # which is max
        j_array = np.where(j_array== major, 1, 0) # if is major then
        
        
        # generate the frequencies of all 4 alleles
        ld_ops = np.zeros(shape=(4),dtype=(np.float32)) # as weighting is fractional
        for k in range(0,nSeqs): # for each sequence, see which bin it fits in. then rather than inc by 1 . increment by weighting
            
        








        
