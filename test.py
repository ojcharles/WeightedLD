import unittest
import io                   # for capturing stdout prints from ld
import sys
import WeightedLD as wld    # import the program as module
import numpy as np
import re


class TestStuff(unittest.TestCase):
    min_acgt = 0.8
    min_variability = 0.02

    def test_read_fasta(self):
        # test the fasta read and hashmap to integers
        file = "tests/t1_henikoff_paper.fasta"
        alignment = wld.read_fasta(file)
        self.assertEqual(alignment.sum(), 65)

    def test_var_sitesHK(self):
        file = "tests/t1_henikoff_paper.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        # 0: >80% ambig ->False, 1: >80% indel -> False
        self.assertEqual(var_sites_HK.tolist(), [
                         False, False, True, True, True, True, True])

    def test_var_sitesLD(self):
        # var_stes_ld should have 1 extra False
        file = "tests/t6_varsites_hk_ld.fasta"
        alignment = wld.read_fasta(file)
        min_variability = 0.2  # high
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, min_variability)
        self.assertNotEqual(var_sites_HK[1], var_sites_LD[1])

    def test_hkw_simple(self):
        # should return weights as in the 1994 paper, but normal. ignores first two sites as in [test_var_sites]
        file = "tests/t1_henikoff_paper.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        alignment = alignment[:, var_sites_HK]
        weightsHK = wld.henikoff_weighting(alignment)
        self.assertTrue(np.allclose(weightsHK,
                                    np.array([0.5, 0.5, 0.5, 0.5, 1.0]),
                                    rtol=1e-02, atol=1e-02))

    def test_hkw_complex(self):
        # seq0 has a unique base at site1, so is most distinct, ensure this
        file = "tests/t2_henikoff_complex1.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        alignment = alignment[:, var_sites_HK]
        weightsHK = wld.henikoff_weighting(alignment)
        self.assertEqual(weightsHK[0], 1.0)

    def test_hkw_complex_indel(self):
        # checks indels are handled correctly, here last seq is most unique and contains two insertions
        file = "tests/t3_henikoff_complex2.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        alignment = alignment[:, var_sites_HK]
        weightsHK = wld.henikoff_weighting(alignment)
        self.assertEqual(weightsHK[7], 1.0)

    def test_0ld_flatw(self):
        # flat weights, no LD, simple example
        file = "tests/t4_weights1_ld0.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, 0.99, self.min_variability)  # ignore indel
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map)  # run function as normal
        sys.stdout = sys.__stdout__             # Reset redirect.
        self.assertEqual(capturedOutput.getvalue()[
            22:25], "0.0")     # D

    def test_wld_flatw(self):
        # flat weights, no LD, simple example - but one sequence is now weighted -> not 0 LD
        file = "tests/t4_weights1_ld0.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, 0.1, 0.2)  # ignore indel in LD
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map)  # run function as normal
        sys.stdout = sys.__stdout__             # Reset redirect.
        self.assertNotEqual(capturedOutput.getvalue()[
            22:25], "0.0")     # D

    def test_ld_flatw(self):
        # flat weights, total LD, simple example
        file = "tests/t5_weights1_ld0.25.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map)  # run function as normal
        sys.stdout = sys.__stdout__             # Reset redirect.

        self.assertEqual(capturedOutput.getvalue()[22:27], "-0.25")     # D
        self.assertEqual(capturedOutput.getvalue()[32:33], "1")         # r2


""" to write
    def test_ld_weights_complex(self):
        # flat weights, total LD, simple example
        min_acgt = 0.8
        min_variability = 0.02
        file = "tests/t3_henikoff_complex2.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, min_acgt, min_variability)
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]

        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map)  # run function as normal
        sys.stdout = sys.__stdout__             # Reset redirect.

        # handle character output to something structured
        outTable = re.split(r'\n+', capturedOutput.getvalue())
        outTable = [x.split("\t") for x in outTable]
        t = [el[4] for el in outTable[1:7]]
        t = np.array(t, dtype=float)
#        self.assertTrue(np.allclose(np.array(t),
#                                    np.array(
#                                        ['0.0915', '0.0824', '0.0824', '0.1586', '0.1586', '1.0'], dtype=float)
#                                    ))
"""


def test_vcf(self):
    # this code tests the whole pipeline with a vcf file, as in the 1000 genomes vcf v4.2
    filename = "tests/t7_1000genome.vcf"
    alignment, site_map = wld.handle_vcf(filename)
    weights = wld.henikoff_weighting(alignment)
    wld.ld(alignment, weights, site_map)

    self.assertEqual(round(weights.mean(), 3), 0.002)


if __name__ == '__main__':
    unittest.main()
