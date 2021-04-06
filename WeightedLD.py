# Oscar Charles 210404
# #Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 
# using matt hahn molpopgen book
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
        
        #----- remove any sequences with indel character
        # which index values need removing?
        i_array = a.alignment_array[:,i]
        remove = np.where(i_array == 4) # which seqs contain illegal char in i
        j_array = a.alignment_array[:,j]
        remove = np.append(remove,np.where(j_array == 4)) # which sequences in total contain illegals?
        remove = np.unique(remove)
        
        # remove any sequence as identified above, which contained a illegal in i or j
        i_array = np.delete(i_array, remove)
        j_array = np.delete(j_array, remove)
        tSeqs = len(i_array) # for this comparison, how many seqs are there?
        
        
        
        
        
        # first site
        # find major allele
        unique_elements, counts_elements = np.unique(i_array, return_counts=True)
        major = unique_elements[counts_elements.argmax()] # which value is max
        i_array = np.where(i_array == major, 1, 0) # if is major value then 1 else 0
        PA = np.count_nonzero(i_array == 1) / tSeqs
        Pa = np.count_nonzero(i_array == 0) / tSeqs
        
        
        
        
        
        # second site
        # find major allelle
        unique_elements, counts_elements = np.unique(j_array, return_counts=True)
        major = unique_elements[counts_elements.argmax()] # which is max
        j_array = np.where(j_array == major, 1, 0) # if is major then
        PB = np.count_nonzero(j_array == 1) / tSeqs
        Pb = np.count_nonzero(j_array == 0) / tSeqs
        
        
        # predicted allele freqs if in equilibrium
        PAB = PA * PB
        PAb = PA * Pb
        PaB = Pa * PB
        Pab = Pa * Pb

        
        
        
        # observed allele freqs
        weights = np.zeros(shape=(tSeqs))
        weights[weights == 0] = 1
        ld_ops = np.zeros(shape=(4),dtype=(np.float32)) # as weighting is fractional
        for k in range(0,tSeqs): # for each sequence, see which bin it fits in. then rather than inc by 1 . increment by weighting
            if i_array[k] == 0 and j_array[k] == 0:
                ld_ops[0] = ld_ops[0] + weights[k]
            elif i_array[k] == 1 and j_array[k] == 1:
                ld_ops[3] = ld_ops[3] + weights[k]
            elif i_array[k] == 0 and j_array[k] == 1:
                ld_ops[1] = ld_ops[1] + weights[k]
            elif i_array[k] == 1 and j_array[k] == 0:
                ld_ops[2] = ld_ops[2] + weights[k]
            else:
                print(k)
            ld_ops = ld_ops / tSeqs
                
        # the vector is now as in the hahn molpopgen book and can be used to generate D values
        # ld_ops is pAB, p
        tD = np.zeros(shape=(4),dtype=(np.float32))
        tD[0] = abs(ld_ops[3] - PAB) # MajMaj
        tD[1] = abs(ld_ops[0] - Pab) # minmin  
        tD[2] = abs(ld_ops[2] - PAb) #Majmin
        tD[3] = abs(ld_ops[1] - PaB) #minMaj
        
        out = (tD[0] + tD[1] + tD[2] + tD[3]) / 4
        out = round(out, 6)
        return out
        
        








        
