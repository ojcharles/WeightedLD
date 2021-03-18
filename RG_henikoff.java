import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.Nucleotides;
import pal.datatype.SpecificAminoAcids;

/**
*
* @author rgoldst
*/
public class HenikofWeighting {
    
    int nSites;
    int nSequences;
    HashMap<Character, Integer> convertBaseToInt = new HashMap<>();
    String[] names;
    char[][] sequences = null;
    
    double minACGT = 0.2;          // Minimum fractions of ACTG at a given site for the site to be included in calculation
    // increase this to remove more noise say 0.5
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        HenikofWeighting heni = new HenikofWeighting();
        heni.run(args[0]); // alignment
    }
    
    
    HenikofWeighting() { // data type hashmap
        convertBaseToInt.put('A', 0);
        convertBaseToInt.put('C', 1);
        convertBaseToInt.put('G', 2);
        convertBaseToInt.put('T', 3);
        convertBaseToInt.put('a', 0);
        convertBaseToInt.put('c', 1);
        convertBaseToInt.put('g', 2);
        convertBaseToInt.put('t', 3);
        convertBaseToInt.put('-', 4);
        convertBaseToInt.put('.', 4);
    }
    
    void run(String fileName) {
        readAlignment(fileName); // read alignment
        henikoff();
    }

    void readAlignment(String fileName) {
        Alignment align = null;
        try {
            Nucleotides dataType = new Nucleotides();
            FileReader in = new FileReader(fileName);
            align = AlignmentReaders.readFastaSequences(in, dataType);
            in.close();
        } catch (IOException e) {
            System.out.println("Error: File not found (IO error)");
            System.exit(1);
        }
        nSites = align.getSiteCount();
        nSequences = align.getSequenceCount();
        names = new String[nSequences];
        System.out.println("Number of sites: " + nSites + "\tNumber of sequences: " + nSequences);
        sequences = new char[nSequences][nSites]; // character matrix
        for (int iSeq = 0; iSeq < nSequences; iSeq++) {
            names[iSeq] = align.getIdentifier(iSeq).toString();
            sequences[iSeq] = align.getAlignedSequenceString(iSeq).toCharArray();
        }
    }

    
    void henikoff() {
        double[] weights = new double[nSequences];      // Resulting sequence weights
        double[] fracACTG = new double[nSequences];     // Avg frac of ACTG in counted sites
        double[] fracOK = new double[nSequences];       // Avg frac of ACTG- in counted istes
        double nSitesCounted = 0.;                      // Number of counted sites
        
        for (int iSite = 0; iSite < nSites; iSite++) {  // Look over all sites
            
            int[] iBaseBySeq = new int[nSequences];     // What base does each sequence have?
            boolean[] okBase = new boolean[nSequences]; // Is the base a canonical base or gap?
            double[] countBase = new double[5];         // How many of each type of base
            double nOK = 0.0;                           // How many sequences have a canonical base or gap
            
            Arrays.fill(iBaseBySeq, 4);                 // Assume gap
            Arrays.fill(okBase, false);                 // Assume not canonical base
            
            // Loop over sequences filling in iBasePerSeq and countBase
            for (int iSeq = 0; iSeq < nSequences; iSeq++) {
                if (convertBaseToInt.containsKey(sequences[iSeq][iSite])) {
                    int iBase = convertBaseToInt.get(sequences[iSeq][iSite]);
                    okBase[iSeq] = true;
                    nOK++;
                    iBaseBySeq[iSeq] = iBase;
                    countBase[iBase]++;
                }
            }
            
            // Compute nBases, the number of different bases at that site
            double nBases = 0.;
            if (countBase[0] > 0.1) nBases++;
            if (countBase[1] > 0.1) nBases++;
            if (countBase[2] > 0.1) nBases++;
            if (countBase[3] > 0.1) nBases++;
            if (countBase[4] > 0.1) nBases++;
            double nRealBase = countBase[0] + countBase[1] + countBase[2] + countBase[3];
            
            // Update weights with the calculation for this site
            if (nBases > 0.1 && (nRealBase/nSequences >= minACGT)) { // Only do this calculation if nBases > 0 and if the fraction of {A,G,C,T} is >= minACGT;
                nSitesCounted++;
                double[] avgWeight = new double[2];
                for (int iSeq = 0; iSeq < nSequences; iSeq++) {
                    if (okBase[iSeq]) { // is actg or -
                        double siteContribution = 1.0/(nBases * countBase[iBaseBySeq[iSeq]]); // contribution to the weight from this site, this is the henikoff
                        weights[iSeq] += siteContribution;
                        avgWeight[0] += siteContribution;
                        avgWeight[1] ++;
                        fracOK[iSeq]++;
                        if (iBaseBySeq[iSeq] < 4) {
                            fracACTG[iSeq]++;
                        }
                    }
                }
                avgWeight[0] /= avgWeight[1];
                for (int iSeq = 0; iSeq < nSequences; iSeq++) {
                    if (!okBase[iSeq]) { // not okbase
                        weights[iSeq] += avgWeight[0]; // don't downweight becuase the sequence was poor quality/ ambigipus add the mean weight.
                        // i.e. the sequence weight is unaffected by ambigious sites. 
                    }
                }                
                
            }   // End of loop over bases
            
        }   // End of loop over sites
        
        // Normalise weights and output
        double norm = 0.0;
        for (int iSeq = 0; iSeq < nSequences; iSeq++) {
            norm = Math.max(norm, weights[iSeq]); // Take the largest weight normalise everything to 1 relative to this.
        }

        for (int iSeq = 0; iSeq < nSequences; iSeq++) { // print output
            System.out.format("%s\t%.4f\t%.4f\t%.4f\n", names[iSeq].replaceAll(" ","_"), weights[iSeq]/norm,
                   (fracACTG[iSeq]/nSitesCounted), (fracOK[iSeq]/nSitesCounted));
        }
        
    }

    
}
