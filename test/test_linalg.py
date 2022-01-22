"""Tests for distos.linalg
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.linalg import dot_product, div
from dictos.error.linear_algebra import InconsistentDataSetError
from gen import random_int, random_string


class LinalgTest(unittest.TestCase):
    def test_linalg_div(self):
        """test suite for linalg.dot_product"""

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            numer = random_int(-n, n)
            denom = random.randint(1, 100)
            with self.subTest(
                "dev of %d elements numerator and %d" % (len(numer), denom)
            ):
                expected = sum([n for n in numer]) / denom
                actural = div(sum([n for n in numer]), denom)
                self.assertEqual(expected, actural)

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            numer = [sp.symbols(random_string(2)) for _ in range(num_max)]
            denom = sp.symbols(random_string(2))
            with self.subTest(
                "dev of %d elements numerator and %s" % (len(numer), denom)
            ):
                expected = sum([n for n in numer]) / denom
                actural = div(sum([n for n in numer]), denom)
                self.assertEqual(expected, actural)

    def test_linalg_exception(self):
        """test suite for exception in linalg"""
        numer = [sp.core.Number(i) for i in range(4)]
        f_set = sp.symbols("f_0:{:d}".format(len(numer) + 1))
        with self.subTest("dot_product with inconsistent numer and f_set"):
            with self.assertRaises(InconsistentDataSetError):
                dot_product(numer, f_set)


if __name__ == "__main__":
    unittest.main()
