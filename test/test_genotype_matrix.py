"""
Test cases for concatenation stage of EVCouplings complex pipeline
Author:
    Anna G. Green
"""

import unittest
import pandas as pd
import numpy as np
import sys
sys.path.append("/Users/agreen/codebases/GenotypeTools")
from genotype_matrix import *
from unittest import TestCase

class TestGenotypeMatrix(TestCase):

    def __init__(self, *args, **kwargs):
        super(TestGenotypeMatrix, self).__init__(*args, **kwargs)

        self.matrix_mapped = np.array([
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 4],
            [1, 2, 3, 4, 0]
        ])

        self.matrix = np.array([
            ["-", "A", "C", "G", "T"],
            ["-", "A", "C", "G", "T"],
            ["-", "A", "C", "G", "T"],
            ["A", "C", "G", "T", "-"]
        ])

        self.positions = ['123', '124', '150', '176', '189']

        self.strains = ["first", "second", "third", "fourth"]

        self.matrix_file= "genotype_matrix.csv"
        self.genotype_data = pd.read_csv("genotype_data.csv")
        self.frequencies = np.array([
            [0.75, 0.25, 0.0, 0.0, 0.0],
            [0.0, 0.75, 0.25, 0.0, 0.0],
            [0.0, 0.0, 0.75, 0.25, 0.0],
            [0.0, 0.0, 0.0, 0.75, 0.25],
            [0.25, 0.0, 0.0, 0.0, 0.75],
        ])

    def test_calculate_frequencies(self):
        frequencies = calculate_frequencies(self.matrix_mapped, 5)
        self.assertEqual(frequencies.tolist(), self.frequencies.tolist())

    def test_from_df(self):
        gm = GenotypeMatrix.from_df(self.matrix_file)
        self.assertEqual(self.strains, gm.strains)
        self.assertEqual(self.positions, gm.positions)
        self.assertEqual(self.matrix.tolist(), gm.matrix.tolist())

    def test_make_matrix(self):
        gm=GenotypeMatrix(self.strains,[int(x) for x in self.positions],genotype_data=self.genotype_data)
        gm.make_matrix()
        print(gm.matrix, gm.strains, gm.positions)
        self.assertEqual(gm.matrix.tolist(), self.matrix.tolist())

    def test_major_alleles(self):
        gm = GenotypeMatrix.from_df(self.matrix_file)
        print(gm.major_alleles)
        self.assertEqual([0,1,2,3,4],gm.major_alleles.tolist())

    def test_write_fasta(self):
        pass

    def test_write_GEMMA(self):
        pass

    def test_to_df(self):
        pass
if __name__ == '__main__':
    unittest.main()
