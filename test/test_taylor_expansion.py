"""Tests for distos.lagrangian_polynomial
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import numpy as np
import random

from dictos.taylor_expansion import taylor_series, derivative_symbol


class UtilsTest(unittest.TestCase):
    def test_derivative_symbol(self):
        """
        test suite for taylor_expansion.derivative_symbol.
        """

        for _ in range(10):
            n = random.randint(0, np.iinfo(np.int32).max)
            with self.subTest(n):
                expected = sp.symbols("f^(%d)" % n)
                acctual = derivative_symbol("f", n)
                self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
