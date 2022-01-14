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
                actual = lagrangian_basis(x, degree=1, point_at=i)
                self.assertEqual(expected[i], actual)

        expected = [
            (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2)),
            (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2)),
            (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1)),
        ]
        for i in range(3):
            with self.subTest("3-point formulation  defined at x_set[%d]" % i):
                actual = lagrangian_basis(x, degree=2, point_at=i)
                self.assertEqual(expected[i], actual)

        expected = [
            (x - x1) * (x - x2) * (x - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)),
            (x - x0) * (x - x2) * (x - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)),
            (x - x0) * (x - x1) * (x - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)),
            (x - x0) * (x - x1) * (x - x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2)),
        ]
        for i in range(4):
            with self.subTest("4-point formulation defined at x_set[%d]" % i):
                actual = lagrangian_basis(x, degree=3, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=4, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=5, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=6, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=7, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=8, point_at=i)
                self.assertEqual(expected[i], actual)

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
                actual = lagrangian_basis(x, degree=9, point_at=i)
                self.assertEqual(expected[i], actual)

    def test_lagrangian_poly(self):
        """test suite for lagrangian_polynomial.lagrangian_poly."""
        pass
        x = sp.symbols("x")
        h = sp.symbols("h")

        f0 = sp.symbols("f0")
        f1 = sp.symbols("f1")
        f2 = sp.symbols("f2")
        f3 = sp.symbols("f3")
        f4 = sp.symbols("f4")

        # subtests for forward formulation
        expected = [
            -f0 * (-h + x) / h + f1 * x / h,
            f0 * (-2 * h + x) * (-h + x) / (2 * h ** 2)
            - f1 * x * (-2 * h + x) / h ** 2
            + f2 * x * (-h + x) / (2 * h ** 2),
            -f0 * (-3 * h + x) * (-2 * h + x) * (-h + x) / (6 * h ** 3)
            + f1 * x * (-3 * h + x) * (-2 * h + x) / (2 * h ** 3)
            - f2 * x * (-3 * h + x) * (-h + x) / (2 * h ** 3)
            + f3 * x * (-2 * h + x) * (-h + x) / (6 * h ** 3),
            f0 * (-4 * h + x) * (-3 * h + x) * (-2 * h + x) * (-h + x) / (24 * h ** 4)
            - f1 * x * (-4 * h + x) * (-3 * h + x) * (-2 * h + x) / (6 * h ** 4)
            + f2 * x * (-4 * h + x) * (-3 * h + x) * (-h + x) / (4 * h ** 4)
            - f3 * x * (-4 * h + x) * (-2 * h + x) * (-h + x) / (6 * h ** 4)
            + f4 * x * (-3 * h + x) * (-2 * h + x) * (-h + x) / (24 * h ** 4),
        ]
        point = [2, 3, 4, 5]
        for i in range(len(expected)):
            with self.subTest("%d-point forward formulation" % point[i]):
                x_set = [j * h for j in range(point[i])]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = lagrangian_poly(x, x_set, f_set)

                self.assertEqual(expected[i], actual)

        # subtests for backward formulation
        expected = [
            -f0 * x / h + f1 * (h + x) / h,
            f0 * x * (h + x) / (2 * h ** 2)
            - f1 * x * (2 * h + x) / h ** 2
            + f2 * (h + x) * (2 * h + x) / (2 * h ** 2),
            -f0 * x * (h + x) * (2 * h + x) / (6 * h ** 3)
            + f1 * x * (h + x) * (3 * h + x) / (2 * h ** 3)
            - f2 * x * (2 * h + x) * (3 * h + x) / (2 * h ** 3)
            + f3 * (h + x) * (2 * h + x) * (3 * h + x) / (6 * h ** 3),
            f0 * x * (h + x) * (2 * h + x) * (3 * h + x) / (24 * h ** 4)
            - f1 * x * (h + x) * (2 * h + x) * (4 * h + x) / (6 * h ** 4)
            + f2 * x * (h + x) * (3 * h + x) * (4 * h + x) / (4 * h ** 4)
            - f3 * x * (2 * h + x) * (3 * h + x) * (4 * h + x) / (6 * h ** 4)
            + f4 * (h + x) * (2 * h + x) * (3 * h + x) * (4 * h + x) / (24 * h ** 4),
        ]
        point = [2, 3, 4, 5]
        for i in range(len(expected)):
            with self.subTest("%d-point backward formulation" % point[i]):
                x_set = [j * h for j in range(-point[i] + 1, 1)]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = lagrangian_poly(x, x_set, f_set)

                self.assertEqual(expected[i], actual)

        # subtests for central formulation
        expected = [
            f0 * x * (-h + x) / (2 * h ** 2)
            - f1 * (-h + x) * (h + x) / h ** 2
            + f2 * x * (h + x) / (2 * h ** 2),
            f0 * x * (-2 * h + x) * (-h + x) * (h + x) / (24 * h ** 4)
            - f1 * x * (-2 * h + x) * (-h + x) * (2 * h + x) / (6 * h ** 4)
            + f2 * (-2 * h + x) * (-h + x) * (h + x) * (2 * h + x) / (4 * h ** 4)
            - f3 * x * (-2 * h + x) * (h + x) * (2 * h + x) / (6 * h ** 4)
            + f4 * x * (-h + x) * (h + x) * (2 * h + x) / (24 * h ** 4),
        ]
        point = [3, 5]
        for i in range(len(expected)):
            with self.subTest("%d-point forward formulation" % point[i]):
                begin_ = -(point[i] - 1) // 2
                end_ = (point[i] - 1) // 2 + 1
                stencil = range(begin_, end_)
                x_set = [j * h for j in stencil]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = lagrangian_poly(x, x_set, f_set)

                self.assertEqual(expected[i], actual)

    def test_derivative(self):
        """test suite for lagrangian.polynomial.derivative."""
        x = sp.symbols("x")

        # nth derivative of sin(x)
        expr = sp.sin(x)
        expected = [1, 0, -1, 0, 1]
        for i in range(len(expected)):
            n = i + 1
            with self.subTest("%d-th derivative of sin(x) at x=0" % n):
                actual = derivative(expr, x, deriv=n)

                self.assertEqual(expected[i], actual)

        # nth derivative of exp(x)
        expr = sp.exp(x)
        for i in range(5):
            n = i + 1
            with self.subTest("%d-th derivative of exp(x) at x=0" % n):
                actual = derivative(expr, x, deriv=n)

                self.assertEqual(1, actual)

        # subtests for finite difference
        h = sp.symbols("h")
        f0 = sp.symbols("f0")
        f1 = sp.symbols("f1")
        f2 = sp.symbols("f2")
        f3 = sp.symbols("f3")
        f4 = sp.symbols("f4")

        # subtests for forward difference
        expected = [
            (-f0 + f1) / h,
            (-3 * f0 + 4 * f1 - f2) / (2 * h),
            (-11 * f0 + 18 * f1 - 9 * f2 + 2 * f3) / (6 * h),
            (-25 * f0 + 48 * f1 - 36 * f2 + 16 * f3 - 3 * f4) / (12 * h),
        ]
        point = [2, 3, 4, 5]
        for i in range(len(expected)):
            with self.subTest(
                "%d-point forward difference for 1st derivative" % point[i]
            ):
                x_set = [j * h for j in range(point[i])]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = derivative(lagrangian_poly(x, x_set, f_set), x, deriv=1)

                self.assertEqual(expected[i], actual)

        # subtests for backward difference
        expected = [
            (-f0 + f1) / h,
            (f0 - 4 * f1 + 3 * f2) / (2 * h),
            (-2 * f0 + 9 * f1 - 18 * f2 + 11 * f3) / (6 * h),
            (3 * f0 - 16 * f1 + 36 * f2 - 48 * f3 + 25 * f4) / (12 * h),
        ]
        point = [2, 3, 4, 5]
        for i in range(len(expected)):
            with self.subTest(
                "%d-point backward difference for 1st derivative" % point[i]
            ):
                x_set = [j * h for j in range(-point[i] + 1, 1)]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = derivative(lagrangian_poly(x, x_set, f_set), x, deriv=1)

                self.assertEqual(expected[i], actual)

        # subtests for central difference
        expected = [
            (-f0 + f2) / (2 * h),
            (f0 - 8 * f1 + 8 * f3 - f4) / (12 * h),
        ]
        point = [3, 5]
        for i in range(len(expected)):
            with self.subTest(
                "%d-point forward difference for 1st derivative" % point[i]
            ):
                begin_ = -(point[i] - 1) // 2
                end_ = (point[i] - 1) // 2 + 1
                stencil = range(begin_, end_)
                x_set = [j * h for j in stencil]
                f_set = sp.symbols("f0:{:d}".format(len(x_set)))
                actual = derivative(lagrangian_poly(x, x_set, f_set), x, 1)

                self.assertEqual(expected[i], actual)


if __name__ == "__main__":
    unittest.main()
