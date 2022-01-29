"""Tests for distos.error.taylor_expansion
"""
import sys

sys.path.insert(1, "..")

import unittest

from dictos.taylor_expansion import taylor_series, derivative_symbol
from dictos.error.taylor_expansion import (
    UnsupportedOrderOfDerivativeError,
    NumberOfExpansionTermsIsNotNaturalNumberError,
)


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

    def test_error_taylor_expansion_exception(self):
        """
        test suite for exceptions in taylor_expansion.
        """

        up_to = -1
        with self.subTest(f"taylor_series with invalid up_to {up_to}"):
            with self.assertRaises(NumberOfExpansionTermsIsNotNaturalNumberError):
                taylor_series("h", up_to)

        deriv = -1
        with self.subTest(f"symbols with invalid deriv {deriv}"):
            with self.assertRaises(UnsupportedOrderOfDerivativeError):
                derivative_symbol("h", deriv)


if __name__ == "__main__":
    unittest.main()
