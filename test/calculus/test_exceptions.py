"""Tests for distos.error.finite_difference
"""

import sys

sys.path.insert(1, "..")

import unittest

from dictos.calculus.finite_difference import coefficients
from dictos.calculus.exceptions import UnsupportedOrderOfDerivativeError


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

    def test_error_finite_difference_coefficients_exception(self):
        """
        test suite for finite_difference.coefficients exceptions.
        """

        for deriv in range(-10, 1):
            with self.subTest(f"{deriv}-th order of derivative"):
                with self.assertRaises(UnsupportedOrderOfDerivativeError):
                    coefficients([-2, -1, 0, 1, 2], deriv)


if __name__ == "__main__":
    unittest.main()
