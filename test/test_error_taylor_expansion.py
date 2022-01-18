"""Tests for distos.error.taylor_expansion
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.error.taylor_expansion import UnsupportedOrderOfDerivativeError


class ErrorTaylorExpansionTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_taylor_expansion_UnsupportedOrderOfDerivativeError(self):
        """
        test suite for error.taylor_expansion.UnsupportedOrderOfDerivativeError.
        """

        deriv = -1
        with self.subTest(deriv):
            if deriv < 0:
                raise UnsupportedOrderOfDerivativeError(deriv)


if __name__ == "__main__":
    unittest.main()
