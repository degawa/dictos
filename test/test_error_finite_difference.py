"""Tests for distos.error.finite_difference
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.error.finite_difference import UnsupportedOrderOfDerivativeError


class ErrorFiniteDifferenceTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_finite_difference_UnsupportedOrderOfDerivativeError(self):
        """
        test suite for error.finite_difference.UnsupportedOrderOfDerivativeError.
        """

        deriv = 0
        with self.subTest(deriv):
            if deriv < 1:
                raise UnsupportedOrderOfDerivativeError(deriv)


if __name__ == "__main__":
    unittest.main()
