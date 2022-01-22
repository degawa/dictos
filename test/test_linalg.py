"""Tests for distos.linalg
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.linalg import dot_product
from dictos.error.linear_algebra import InconsistentDataSetError


class LinalgTest(unittest.TestCase):
    def test_linalg_exception(self):
        """test suite for exception in linalg"""
        numer = [sp.core.Number(i) for i in range(4)]
        f_set = sp.symbols("f_0:{:d}".format(len(numer) + 1))
        with self.subTest("dot_product with inconsistent numer and f_set"):
            with self.assertRaises(InconsistentDataSetError):
                dot_product(numer, f_set)


if __name__ == "__main__":
    unittest.main()
