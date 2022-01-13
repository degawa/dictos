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

        with self.subTest(1):
            expected = (x - x1) / (x0 - x1)
            acctual = lagrangian_basis(x, degree=1, point_at=0)
            self.assertEqual(expected, acctual)

            expected = (x - x0) / (x1 - x0)
            acctual = lagrangian_basis(x, degree=1, point_at=1)
            self.assertEqual(expected, acctual)

        with self.subTest(2):
            expected = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))
            acctual = lagrangian_basis(x, degree=2, point_at=0)
            self.assertEqual(expected, acctual)

            expected = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))
            acctual = lagrangian_basis(x, degree=2, point_at=1)
            self.assertEqual(expected, acctual)

            expected = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1))
            acctual = lagrangian_basis(x, degree=2, point_at=2)
            self.assertEqual(expected, acctual)

        with self.subTest(3):
            expected = (
                (x - x1) * (x - x2) * (x - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3))
            )
            acctual = lagrangian_basis(x, degree=3, point_at=0)
            self.assertEqual(expected, acctual)

            expected = (
                (x - x0) * (x - x2) * (x - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3))
            )
            acctual = lagrangian_basis(x, degree=3, point_at=1)
            self.assertEqual(expected, acctual)

            expected = (
                (x - x0) * (x - x1) * (x - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3))
            )
            acctual = lagrangian_basis(x, degree=3, point_at=2)
            self.assertEqual(expected, acctual)

        with self.subTest(4):
            expected = (
                (x - x1)
                * (x - x2)
                * (x - x3)
                * (x - x4)
                / ((x0 - x1) * (x0 - x2) * (x0 - x3) * (x0 - x4))
            )
            acctual = lagrangian_basis(x, degree=4, point_at=0)
            self.assertEqual(expected, acctual)

            expected = (
                (x - x0)
                * (x - x2)
                * (x - x3)
                * (x - x4)
                / ((x1 - x0) * (x1 - x2) * (x1 - x3) * (x1 - x4))
            )
            acctual = lagrangian_basis(x, degree=4, point_at=1)
            self.assertEqual(expected, acctual)

            expected = (
                (x - x0)
                * (x - x1)
                * (x - x3)
                * (x - x4)
                / ((x2 - x0) * (x2 - x1) * (x2 - x3) * (x2 - x4))
            )
            acctual = lagrangian_basis(x, degree=4, point_at=2)
            self.assertEqual(expected, acctual)

            expected = (
                (x - x0)
                * (x - x1)
                * (x - x2)
                * (x - x4)
                / ((x3 - x0) * (x3 - x1) * (x3 - x2) * (x3 - x4))
            )
            acctual = lagrangian_basis(x, degree=4, point_at=3)
            self.assertEqual(expected, acctual)

        with self.subTest(5):
            expected = (
                (x - x1)
                * (x - x2)
                * (x - x3)
                * (x - x4)
                * (x - x5)
                / ((x0 - x1) * (x0 - x2) * (x0 - x3) * (x0 - x4) * (x0 - x5))
            )
            acctual = lagrangian_basis(x, degree=5, point_at=0)
            self.assertEqual(expected, acctual)

        with self.subTest(6):
            expected = (
                (x - x1)
                * (x - x2)
                * (x - x3)
                * (x - x4)
                * (x - x5)
                * (x - x6)
                / (
                    (x0 - x1)
                    * (x0 - x2)
                    * (x0 - x3)
                    * (x0 - x4)
                    * (x0 - x5)
                    * (x0 - x6)
                )
            )
            acctual = lagrangian_basis(x, degree=6, point_at=0)
            self.assertEqual(expected, acctual)

        with self.subTest(7):
            expected = (
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
            )
            acctual = lagrangian_basis(x, degree=7, point_at=0)
            self.assertEqual(expected, acctual)

        with self.subTest(8):
            expected = (
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
            )
            acctual = lagrangian_basis(x, degree=8, point_at=0)
            self.assertEqual(expected, acctual)

        with self.subTest(9):
            expected = (
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
            )
            acctual = lagrangian_basis(x, degree=9, point_at=0)
            self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
