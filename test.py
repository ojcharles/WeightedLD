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

        self.assertEqual(alignment.sum(), 76)

    def test_var_sitesHK(self):
        file = "tests/t1_henikoff_paper.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)

        # non-variable site
        self.assertEqual(var_sites_HK[0], False)
        # too many non atgc
        self.assertEqual(var_sites_HK[1], False)
        # 0.8 are not acgt and we use > not >=
        self.assertEqual(var_sites_HK[2], False)

    def test_var_sitesLD(self):
        # var_stes_ld should have 1 extra False
        file = "tests/t6_varsites_hk_ld.fasta"
        alignment = wld.read_fasta(file)
        min_variability = 0.2
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
        # seq0 has a unique base at site0, so is most distinct, ensure this
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

    def test_0ld_flatw_simple(self):
        # flat weights, no LD, simple example
        file = "tests/t4_weights1_ld0.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]

        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.

        wld.ld(alignment, weightsHK, site_map)  # run function as normal

        sys.stdout = sys.__stdout__             # Reset redirect.
        self.assertEqual(capturedOutput.getvalue()[22:25], "0.0")     # D

    def test_ld_flatw_simple(self):
        # flat weights, total LD, simple example
        file = "tests/t5_weights1_ld0.25.fasta"
        alignment = wld.read_fasta(file)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, self.min_acgt, self.min_variability)
        weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]

        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.

        wld.ld(alignment, weightsHK, site_map)  # run function as normal

        sys.stdout = sys.__stdout__             # Reset redirect.

        self.assertEqual(capturedOutput.getvalue()[22:27], "-0.25")     # D
        self.assertEqual(capturedOutput.getvalue()[32:33], "1")         # r2

    def test_ld_considers_only_Major_donMinor(self):
        # flat weights, total LD, simple example
        min_acgt = 0.8
        min_variability = 0.02
        file = "tests/t7_henikoff_paper_ld.fasta"
        alignment = wld.read_fasta(file)
        # print(alignment)
        var_sites_HK, var_sites_LD = wld.compute_variable_sites(
            alignment, min_acgt, min_variability)
        # print(var_sites_LD)
        # weightsHK = wld.henikoff_weighting(alignment[:, var_sites_HK])
        weightsHK = np.array([1, 1, 1, 1, 1, 1])
        alignment = alignment[:, var_sites_LD]
        site_map = np.where(var_sites_LD)[0]

        # capture stdout
        capturedOutput = io.StringIO()          # Create io object
        sys.stdout = capturedOutput             # to which we redirect stdout.

        wld.ld(alignment, weightsHK, site_map)  # run function as normal

        sys.stdout = sys.__stdout__             # Reset redirect.
        print(capturedOutput.getvalue())


if __name__ == '__main__':
    unittest.main()
