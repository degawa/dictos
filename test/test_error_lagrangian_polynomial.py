"""Tests for distos.error.lagrangian_polynomial
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.error.lagrangian_polynomial import (
    InconsistentDataSetError,
    ViolateDegreeOfPolynomialAssumption,
    InconsistentDataSetAndDegreeOfPolynomialError,
)


class ErrorLagrangianPolynomialTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_stencil_InconsistentDataSetError(self):
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
    def test_error_stencil_ViolateDegreeOfPolynomialAssumption(self):
        """
        test suite for error.lagrangian_polynomial.ViolateDegreeOfPolynomialAssumption.
        """

        degree = 1
        with self.subTest(degree):
            if degree <= 0:
                raise ViolateDegreeOfPolynomialAssumption(degree)

        degree = 0
        with self.subTest(degree):
            if degree <= 0:
                raise ViolateDegreeOfPolynomialAssumption(degree)

    @unittest.expectedFailure
    def test_error_stencil_InconsistentDataSetAndDegreeOfPolynomialError(self):
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


if __name__ == "__main__":
    unittest.main()
