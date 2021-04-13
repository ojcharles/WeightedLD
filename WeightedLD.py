# Oscar Charles 210404
# Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights" 
# using matt hahn molpopgen book
# This code calculates sequence weights using the Henikoff formula from a multiple sequuence alignment
# usage python WeightedLD.py alignment.fasta

#---------- todo
# sequence weighting - treat 4 as a not ok character
# LD - why are some values being returned infinte?

from Bio import AlignIO
import os
import numpy as np
os.chdir('C:\\Oscar\\OneDrive\\UCL\\code\\WeightedLD')

### modifications
minACGT = 0.20   # Minimum fractions of ACTG at a given site for the site to be included in calculation. increase this to remove more noise say 0.5
#msa = sys.argv[0]
#alignmentFile = "all_raw_best_msa_man4.fasta"
#alignmentFile = "test_weights1_LD0.fasta"
alignmentFile = "test_weights1_hahn2.fasta"
### end


class MSA:
    ### this handles the alignment
    # read alignment
    # generate allelle matrix
    # todo
    # check the LD, the weighted values are all weird and wonderful
    def __init__(self, alignmentFile, minACGT):
        self.alignment = AlignIO.read(alignmentFile, "fasta")
        self.minACGT = minACGT
        self.nSeqs = len(self.alignment)
        self.nSites = self.alignment.get_alignment_length()
        self.alignment_array = self.alignment2alignment_array(self.alignment)
        self.var_sites = self.which_pos_var()
        #self.weighting = self.henikoff_weighting(self)
        # LD will be called separately
    
    
    def alignment2alignment_array(self, alignment):
        alignment_array = np.array([list(rec) for rec in alignment], np.dtype("U4")) # alignment as matrix, of dtype unicode
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
           #".":4,- pref behaviour is throw an error
           "k":5,"m":5,"n":5,"r":5,"s":5,"v":5,"y":5,"b":5,"w":5,"h":5,"d":5} # ambiguity codes - impute N or mean score
        # translate character matrix to integer using dict
        alignment_array = np.vectorize(BaseToInt.get)(alignment_array)
        alignment_array = alignment_array.astype(np.uint8()) #convert to 1 byte
        return alignment_array
            
    
    def which_pos_var(self):
        # in
            # alignment array
        # ignores sites where min ACTG fraction not met
        # out
            # vector of variable position ignore character 4
        which_var = np.zeros(shape=(self.nSites),dtype=(np.uint8))
        for pos in range(0,self.nSites):
            array = self.alignment_array[:,pos]
            # if is too little information
            if ( np.count_nonzero(array >= 4) / self.nSeqs) > minACGT: # if too many gaps / ambigious.
                # //todo try something reasonable
                continue # goto next
            a = np.unique(array)
            a = a[a < 4] # remove the indels
            a = len(a)
            if a > 1:
                which_var[pos] = 1
        out = np.where(which_var == 1) # vector of 0 and 1 ->  index of var pos
        return out[0] # force return 1darray
    
    
    def henikoff_weighting(self):
        # todo treat val of 4 as not ok.
        # generate arrays of nSeqs
        # for each pos
            # sum the weighting to seqs
        # then normalise
        weights = np.zeros(shape=(self.nSeqs),dtype=(np.float32()))    # the output
        fracOK = np.zeros(shape=(self.nSeqs),dtype=(np.float32()))     
        nSitesCounted = 0                                           # we might not count all sites in the weighting
        okBaseVals = [0,1,2,3,4]                                      # which values are ok
        
        #---------- loop over each site, and for each seq keep tally of cumulative weighting
        for iSite in self.var_sites:                                   # for each variable site
            #print(iSite)
            array = self.alignment_array[:,iSite]                      # the already converted array i.e. actg--> 01234
            unique_elements, counts_elements = np.unique(array, return_counts=True)
            okBase = np.in1d(array, okBaseVals)                     # vector of T/F for okayness ignores anythin other than actg-.
            #tSeqs = np.count_nonzero(okBase)                        # how many seqs are ok and considered?
            nSitesCounted = nSitesCounted + 1
            countBase = np.zeros(shape=(5),dtype=(np.int16()))
            countBase[0] = np.count_nonzero(array == 0)
            countBase[1] = np.count_nonzero(array == 1)
            countBase[2] = np.count_nonzero(array == 2)
            countBase[3] = np.count_nonzero(array == 3)
            countBase[4] = np.count_nonzero(array == 4) # oscar just assume this is working out for 5 bases not 4.
            tSeqs = countBase[0] + countBase[1] + countBase[2] + countBase[3] + countBase[4]  # //todo same as tSeq
            if ( ( tSeqs - countBase[4] ) / self.nSeqs ) < self.minACGT:   #if too many missing vals then go to next
                continue
            avgWeight = np.zeros(shape=(2),dtype=(np.float32()))
            for iSeq in range(0,self.nSeqs):  # //todo update this so that it defulats to 0 and only checks okbases
                if okBase[iSeq]: # calculate the site contribution
                    iSeq_base = array[iSeq]
                    siteContribution = 1.0/(tSeqs * countBase[iSeq_base]);  # key    contribution to the weight from this site, this is the henikoff
                    weights[iSeq] = weights[iSeq] + siteContribution        # Key    For this seq, keep a tab of its cumulative weight over all sites
                    avgWeight[0] = avgWeight[0] + siteContribution
                    avgWeight[1] = avgWeight[1] + 1
                    fracOK[iSeq] = fracOK[iSeq] + 1
                
            avgWeight[0] =  avgWeight[0] / avgWeight[1]
            # give any imbigious charcters the mean weighting score
            for iSeq in range(0,self.nSeqs):   # if not okbase, give it the average for this pos
                if not okBase[iSeq]: 
                    weights[iSeq] += avgWeight[0]
        # end of loop over sites
        
        
        #---------- normalise
        norm = weights / weights.max()
        fracOK = fracOK / nSitesCounted
        
        #---------- stdout
        print("seq\tfracOK\tWeight")
        for iSeq in range(0,self.nSeqs):
            print(str(iSeq) + "\t" + str(fracOK[iSeq]) + "\t" + str(norm[iSeq]))
            
        #end
        print("alles schon")
        return norm
        
        
    
    
    def LD(self, weights):
        # generate LD metric D and R2
        print("posa\tposb\tD\tD'\tR2") # stdout headers
        iLoops = 0
        outer_loop = self.var_sites[0:len(self.var_sites)-1]
        print(self.var_sites)
        for i in outer_loop:
            iLoops += 1
            inner_loop = self.var_sites[iLoops:len(self.var_sites)]
            for j in inner_loop:
                #print(str(i) + "   " + str(j))
                # compare each pairwise position
                
                #----- remove any sequences with indel or illegal character
                # which index values need removing?
                i_array = self.alignment_array[:,i]
                remove = np.where(i_array >= 4) # which seqs contain illegal char in i
                j_array = self.alignment_array[:,j]
                remove = np.append(remove,np.where(j_array >= 4)) # which sequences in total contain illegals?
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
                ld_ops = np.zeros(shape=(4),dtype=(np.float32)) # as weighting is fractional
                for k in range(0,tSeqs): # for each sequence, see which bin it fits in. then rather than inc by 1 . increment by weighting
                    if i_array[k] == 0 and j_array[k] == 0:     #Oab
                        ld_ops[0] = ld_ops[0] + tWeights[k]
                    elif i_array[k] == 1 and j_array[k] == 1:   #OAB
                        ld_ops[3] = ld_ops[3] + tWeights[k]
                    elif i_array[k] == 0 and j_array[k] == 1:   #OAb
                        ld_ops[1] = ld_ops[1] + tWeights[k]
                    elif i_array[k] == 1 and j_array[k] == 0:   #OaB
                        ld_ops[2] = ld_ops[2] + tWeights[k]
                    else:
                        print(k)
                ld_ops = ld_ops / sum(tWeights)
                
                
                
                
                
                # ----- predicted allele freqs if 0 LD
                PAB = PA * PB
                PAb = PA * Pb
                PaB = Pa * PB
                Pab = Pa * Pb
                
                            
                        
                            
                
                # ----- Caclulate D [Linkage Disequilibrium]
                # the vector is now as in the hahn molpopgen book and can be used to generate D values
                # ld_ops is pAB, p
                tD = np.zeros(shape=(4),dtype=(np.float32))
                tD[0] = abs(ld_ops[3] - PAB) # MajMaj
                tD[1] = abs(ld_ops[0] - Pab) # minmin  
                tD[2] = abs(ld_ops[2] - PAb) #Majmin
                tD[3] = abs(ld_ops[1] - PaB) #minMaj
                
                D = (tD[0] + tD[1] + tD[2] + tD[3]) / 4     #they should be the same anyhow
                
                
                
                # normalised D = D'
                if D < 0:
                    denominator = max([-ld_ops[0],-ld_ops[3]])
                    if denominator == 0:
                        denominator = min([-ld_ops[0],-ld_ops[3]])
                    DPrime = D / denominator  
                else:
                    denominator = min([ld_ops[1],ld_ops[2]])
                    if denominator == 0:
                        denominator = max([ld_ops[1],ld_ops[2]])  
                    DPrime = D / denominator
                    
                
                
                
                # calculate R2
                R2 = D**2 / ( PA * Pa * PB * Pb )
                
                # cat output
                print(str(i)+"\t"+str(j)+"\t"+str(round(D, 4))+"\t"+str(round(DPrime, 4))+"\t"+str(round(R2, 4)))
                


        
        
a = MSA(alignmentFile, minACGT)

weightsHk = a.henikoff_weighting()
weights1 =  np.zeros(shape=(a.nSeqs),dtype=(np.uint16()))
weights1[weights1 == 0] = 1
weights = weights1
ld = a.LD(weightsHk)
    
 
        
