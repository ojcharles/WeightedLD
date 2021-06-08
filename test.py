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
        # example taken from Stephen F. Altschul at NIH NCBI
        file = "tests/t2_henikoff_complex1.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        alignment = alignment[:, var_sites_HK]
        weightsHK = wld.henikoff_weighting(alignment)
        self.assertTrue(np.allclose(weightsHK,
                                    np.array([0.76923077, 0.69230769, 1.0]),
                                    rtol=1e-02, atol=1e-02))

    def test_hkw_complex_indel(self):
        # checks indels are handled correctly, here last seq is most unique and contains two insertions
        file = "tests/t3_henikoff_complex2.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        alignment = alignment[:, var_sites_HK]
        weightsHK = wld.henikoff_weighting(alignment)
        self.assertEqual(weightsHK[7], 1.0)

    def test_0ld_unweighted(self):
        # flat weights, no LD, simple example
        file = "tests/t4_weights1_ld0.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, 0.99, self.min_variability)  # ignore indel
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout, run ld function, split output
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map, r2_threshold=0.0)
        sys.stdout = sys.__stdout__             # Reset redirect.
        out = capturedOutput.getvalue().replace("\n", "\t").split("\t")
        self.assertEqual(out[7], "0.0")     # D
        self.assertEqual(out[8], "0.0")     # DPrime
        self.assertEqual(out[9], "0.0")     # r2

    def test_0ld_weighted(self):
        # flat weights, no LD, simple example - but one sequence is now weighted -> not 0 LD
        file = "tests/t4_weights1_ld0.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, 0.1, 0.2)  # ignore indel in LD
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout, run ld function, split output
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map, r2_threshold=0.0)
        sys.stdout = sys.__stdout__             # Reset redirect.
        out = capturedOutput.getvalue().replace("\n", "\t").split("\t")
        self.assertNotEqual(out[7], "0.0")     # D
        self.assertNotEqual(out[8], "0.0")     # DPrime
        self.assertNotEqual(out[9], "0.0")     # r2

    def test_ld_unweighted(self):
        # flat weights, total LD, simple example
        file = "tests/t5_weights1_ld0.25.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        weightsHK = wld.henikoff_weighting(
            alignment[:, var_sites_HK])  # all equal 1
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout, run ld function, split output
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map, r2_threshold=0.0)
        sys.stdout = sys.__stdout__             # Reset redirect.
        out = capturedOutput.getvalue().replace("\n", "\t").split("\t")
        self.assertEqual(out[7], "-0.25")     # D
        self.assertEqual(out[8], "1.0")     # DPrime
        self.assertEqual(out[9], "1.0")     # r2

    def test_ld_khan_example(self):
        # tests the unweighted LD statistics are correct by a correct worked example
        # https://pbgworks.org/sites/pbgworks.org/files/measuresoflinkagedisequilibrium-111119214123-phpapp01_0.pdf
        file = "tests/t8_ldstats.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        weightsHK = np.ones(2000)  # enforce unit weights
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout, run ld function, split output
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map, r2_threshold=0.0)
        sys.stdout = sys.__stdout__             # Reset redirect.
        out = capturedOutput.getvalue().replace("\n", "\t").split("\t")
        self.assertEqual(out[7], "0.0699")     # D
        self.assertEqual(out[8], "0.4961")     # DPrime
        self.assertEqual(out[9], "0.0924")     # r2

    def test_ld_khan_example(self):
        # A example where we calculated LD with weighting by hand according to our method
        # We used PLINK to check the evenly weighted calculations, then altered the weights
        file = "tests/t9_hand_written.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, 0.001)
        weightsHK = wld.henikoff_weighting(
            alignment[:, var_sites_HK])  # all equal 1
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]
        # capture stdout, run ld function, split output
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.
        wld.ld(alignment, weightsHK, site_map, r2_threshold=0.0)
        sys.stdout = sys.__stdout__             # Reset redirect.
        out = capturedOutput.getvalue().replace("\n", "\t").split("\t")
        self.assertEqual(out[7], "-0.0399")     # D
        self.assertEqual(out[8], "0.1608")     # DPrime
        self.assertEqual(out[9], "0.0259")     # r2


def test_vcf(self):
    # this code tests the whole pipeline with a vcf file, as in the 1000 genomes vcf v4.2
    # this needs improving
    filename = "tests/t7_1000genome.vcf"
    alignment, site_map = wld.handle_vcf(filename)
    weights = wld.henikoff_weighting(alignment)
    wld.ld(alignment, weights, site_map)

    self.assertEqual(round(weights.mean(), 3), 0.002)


if __name__ == '__main__':
    unittest.main()
