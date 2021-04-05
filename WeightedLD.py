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
    alignment = AlignIO.read(alignmentFile, "fasta")
    nSeqs = len(alignment)
    nSites = alignment.get_alignment_length()
    
    def alignment2alignment_array(alignment):
        alignment_array = np.array([list(rec) for rec in alignment], np.character) # alignment as matrix
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
        
    
    
    def LD_D(alignment_array):
        # generate LD metric D
        
        
    


class LD:
    def __init__(self, alignmentFile,title):
        self.title = title
        self.msa = AlignIO.read(alignmentFile, "fasta")
        
        
        

        

a = LD(title = "oscar",alignmentFile = "test.fasta")
print(a.title)


# np play
        
