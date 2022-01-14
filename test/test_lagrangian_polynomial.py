"""Tests for distos.lagrangian_polynomial
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.lagrangian_polynomial import lagrangian_basis, lagrangian_poly, derivative


class UtilsTest(unittest.TestCase):
    def test_lagrangian_basis(self):
        """
        test suite for lagrangian_polynomial.lagrangian_basis.
        """
        x = sp.symbols("x")
        x0 = sp.symbols("x0")
        x1 = sp.symbols("x1")
        x2 = sp.symbols("x2")
        x3 = sp.symbols("x3")
        x4 = sp.symbols("x4")
        x5 = sp.symbols("x5")
        x6 = sp.symbols("x6")
        x7 = sp.symbols("x7")
        x8 = sp.symbols("x8")
        x9 = sp.symbols("x9")

        expected = [(x - x1) / (x0 - x1), (x - x0) / (x1 - x0)]
        for i in range(2):
            with self.subTest("2-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=1, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)),
            (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2)),
            (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)),
        ]
        for i in range(3):
            with self.subTest("3-point formulation  defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=2, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1) * (x - x2) * (x - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)),
            (x - x0) * (x - x2) * (x - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)),
            (x - x0) * (x - x1) * (x - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)),
            (x - x0) * (x - x1) * (x - x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2)),
        ]
        for i in range(4):
            with self.subTest("4-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=3, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            / ((x0 - x1) * (x0 - x2) * (x0 - x3) * (x0 - x4)),
            (x - x0)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            / ((x1 - x0) * (x1 - x2) * (x1 - x3) * (x1 - x4)),
            (x - x0)
            * (x - x1)
            * (x - x3)
            * (x - x4)
            / ((x2 - x0) * (x2 - x1) * (x2 - x3) * (x2 - x4)),
            (x - x0)
            * (x - x1)
            * (x - x2)
            * (x - x4)
            / ((x3 - x0) * (x3 - x1) * (x3 - x2) * (x3 - x4)),
            (x - x0)
            * (x - x1)
            * (x - x2)
            * (x - x3)
            / ((x4 - x0) * (x4 - x1) * (x4 - x2) * (x4 - x3)),
        ]
        for i in range(5):
            with self.subTest("5-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=4, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            * (x - x5)
            / ((x0 - x1) * (x0 - x2) * (x0 - x3) * (x0 - x4) * (x0 - x5))
        ]
        for i in range(1):
            with self.subTest("6-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=5, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            * (x - x5)
            * (x - x6)
            / ((x0 - x1) * (x0 - x2) * (x0 - x3) * (x0 - x4) * (x0 - x5) * (x0 - x6))
        ]
        for i in range(1):
            with self.subTest("7-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=6, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            * (x - x5)
            * (x - x6)
            * (x - x7)
            / (
                (x0 - x1)
                * (x0 - x2)
                * (x0 - x3)
                * (x0 - x4)
                * (x0 - x5)
                * (x0 - x6)
                * (x0 - x7)
            )
        ]
        for i in range(1):
            with self.subTest("8-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=7, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            * (x - x5)
            * (x - x6)
            * (x - x7)
            * (x - x8)
            / (
                (x0 - x1)
                * (x0 - x2)
                * (x0 - x3)
                * (x0 - x4)
                * (x0 - x5)
                * (x0 - x6)
                * (x0 - x7)
                * (x0 - x8)
            )
        ]
        for i in range(1):
            with self.subTest("9-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=8, point_at=i)
                self.assertEqual(expected[i], acctual)

        expected = [
            (x - x1)
            * (x - x2)
            * (x - x3)
            * (x - x4)
            * (x - x5)
            * (x - x6)
            * (x - x7)
            * (x - x8)
            * (x - x9)
            / (
                (x0 - x1)
                * (x0 - x2)
                * (x0 - x3)
                * (x0 - x4)
                * (x0 - x5)
                * (x0 - x6)
                * (x0 - x7)
                * (x0 - x8)
                * (x0 - x9)
            )
        ]
        for i in range(1):
            with self.subTest("10-point formulation defined at x_set[%d]" % i):
                acctual = lagrangian_basis(x, degree=9, point_at=i)
                self.assertEqual(expected[i], acctual)

            )
            acctual = lagrangian_basis(x, degree=9, point_at=0)
            self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
