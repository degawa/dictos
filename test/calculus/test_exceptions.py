"""Tests for distos.calculus.exceptions
"""

import sys

sys.path.insert(1, "..")

import unittest

from dictos.calculus.finite_difference import coefficients, generate
from dictos.calculus.exceptions import (
    UnsupportedOrderOfDerivativeError,
    InvalidOrderOfAccuracyForCentralFormError,
)
from dictos.utilities.spec import is_valid_accuracy_order_for_generating_central_form


class ErrorFiniteDifferenceTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_finite_difference_UnsupportedOrderOfDerivativeError(self):
        """
        test suite for calculus.exceptions.UnsupportedOrderOfDerivativeError.
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

    @unittest.expectedFailure
    def test_error_finite_difference_InvalidOrderOfAccuracyForCentralFormError(self):
        """
        test suite for calculus.exceptions.InvalidOrderOfAccuracyForCentralFormError.
        """

        acc = 3
        with self.subTest(acc):
            if not is_valid_accuracy_order_for_generating_central_form(acc):
                raise InvalidOrderOfAccuracyForCentralFormError(acc)

    def test_error_generate(self):
        """
        test suite for finite_difference.generate exceptions.
        """

        for acc in range(1, 11, 2):
            with self.subTest(
                f"generate {acc}-th order of accuracy central finite difference"
            ):
                with self.assertRaises(InvalidOrderOfAccuracyForCentralFormError):
                    generate(acc=acc, deriv=2)


if __name__ == "__main__":
    unittest.main()
