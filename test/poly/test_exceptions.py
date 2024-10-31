"""Tests for distos.error.lagrangian_polynomial
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.poly.exceptions import (
    InconsistentDataSetError,
    DegreeOfPolynomialIsNotNaturalNumberError,
    InconsistentDataSetAndDegreeOfPolynomialError,
)
from dictos.discrete.exceptions import DuplicatedPointError, TooNarrowError
from dictos.poly.lagrangian_polynomial import lagrangian_basis, lagrangian_poly


class ErrorLagrangianPolynomialTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_lagrangian_polynomial_InconsistentDataSetError(self):
        """
        test suite for error.lagrangian_polynomial.InconsistentDataSetError.
        """

        x = sp.symbols("x")
        x_set = [-x, 0, x]

        f_set = sp.symbols("f0:3")
        with self.subTest(x_set, f_set):
            if len(x_set) != len(f_set):
                raise InconsistentDataSetError(x_set, f_set)

        f_set = sp.symbols("f0:2")
        with self.subTest(x_set, f_set):
            if len(x_set) != len(f_set):
                raise InconsistentDataSetError(x_set, f_set)

    @unittest.expectedFailure
    def test_error_lagrangian_polynomial_ViolateDegreeOfPolynomialAssumption(self):
        """
        test suite for error.lagrangian_polynomial.ViolateDegreeOfPolynomialAssumption.
        """

        degree = 1
        with self.subTest(degree):
            if degree <= 0:
                raise DegreeOfPolynomialIsNotNaturalNumberError(degree)

        degree = 0
        with self.subTest(degree):
            if degree <= 0:
                raise DegreeOfPolynomialIsNotNaturalNumberError(degree)

    @unittest.expectedFailure
    def test_error_lagrangian_polynomial_InconsistentDataSetAndDegreeOfPolynomialError(
        self,
    ):
        """
        test suite for error.lagrangian_polynomial.InconsistentDataSetAndDegreeOfPolynomialError.
        """

        degree = 4
        x_set = [1, 2, 3]
        with self.subTest(degree, x_set):
            if degree != len(x_set) + 1:
                raise InconsistentDataSetAndDegreeOfPolynomialError(degree, x_set)

        degree = 2
        x_set = [1, 2]
        with self.subTest(degree, x_set):
            if degree != len(x_set) + 1:
                raise InconsistentDataSetAndDegreeOfPolynomialError(degree, x_set)

    def test_error_lagrangian_polynomial_exception(
        self,
    ):
        """
        test suite for exception in lagrangian_polynomial.
        """

        x = sp.symbols("x")
        dx = sp.symbols("dx")
        degree = 0
        with self.subTest(f"lagrangian_basis with invalid degree {degree}"):
            with self.assertRaises(DegreeOfPolynomialIsNotNaturalNumberError):
                lagrangian_basis(x, degree, 0)

        degree = 2
        x_set = [0, 0, dx]
        with self.subTest("lagrangian_basis with invalid x_set"):
            with self.assertRaises(DuplicatedPointError):
                lagrangian_basis(x, degree, 0, x_set)

        degree = 3
        x_set = [0, dx, 2 * dx]
        with self.subTest("lagrangian_basis with inconsistent degree and x_set"):
            with self.assertRaises(InconsistentDataSetAndDegreeOfPolynomialError):
                lagrangian_basis(x, degree, 0, x_set)

        degree = 2
        x_set = [0, dx, 2 * dx]
        point_at = -1
        with self.subTest(
            "lagrangian_basis with point_at out of range (less than the lower limit)"
        ):
            with self.assertRaises(ValueError):
                lagrangian_basis(x, degree, point_at, x_set)

        degree = 2
        x_set = [0, dx, 2 * dx]
        point_at = 3
        with self.subTest(
            "lagrangian_basis with point_at out of range (greater than the upper limit)"
        ):
            with self.assertRaises(ValueError):
                lagrangian_basis(x, degree, point_at, x_set)

        x_set = [0, dx, 2 * dx]
        f_set = sp.symbols(f"f0:{(len(x_set) + 1)}")
        with self.subTest("lagrangian_poly with inconsistent x_set and f_set"):
            with self.assertRaises(InconsistentDataSetError):
                lagrangian_poly(x, x_set, f_set)

        x_set = [0]
        f_set = sp.symbols(f"f0:{len(x_set)}")
        with self.subTest("lagrangian_poly with too narrow x_set"):
            with self.assertRaises(TooNarrowError):
                lagrangian_poly(x, x_set, f_set)


if __name__ == "__main__":
    unittest.main()
