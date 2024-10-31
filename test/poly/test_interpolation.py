"""Tests for distos.interplation
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.poly.interpolation import equation, coefficients, truncation_error


class InterpolationTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for interplation.equation.
        """

        f_0 = sp.symbols("f_{0}")
        f_1 = sp.symbols("f_{1}")
        f_2 = sp.symbols("f_{2}")
        f_3 = sp.symbols("f_{3}")
        f_4 = sp.symbols("f_{4}")
        f_5 = sp.symbols("f_{5}")
        f_6 = sp.symbols("f_{6}")
        f_7 = sp.symbols("f_{7}")
        f_8 = sp.symbols("f_{8}")
        f_9 = sp.symbols("f_{9}")
        f_10 = sp.symbols("f_{10}")

        f_m5 = sp.Symbol("f_{-5}")
        f_m4 = sp.Symbol("f_{-4}")
        f_m3 = sp.Symbol("f_{-3}")
        f_m2 = sp.Symbol("f_{-2}")
        f_m1 = sp.Symbol("f_{-1}")
        f_p1 = sp.Symbol("f_{1}")
        f_p2 = sp.Symbol("f_{2}")
        f_p3 = sp.Symbol("f_{3}")
        f_p4 = sp.Symbol("f_{4}")
        f_p5 = sp.Symbol("f_{5}")

        expected = [
            0,
            f_m1 / 2 + f_p1 / 2,
            -f_m2 / 6 + 2 * f_m1 / 3 + 2 * f_p1 / 3 - f_p2 / 6,
            f_m3 / 20
            - 3 * f_m2 / 10
            + 3 * f_m1 / 4
            + 3 * f_p1 / 4
            - 3 * f_p2 / 10
            + f_p3 / 20,
            -f_m4 / 70
            + 4 * f_m3 / 35
            - 2 * f_m2 / 5
            + 4 * f_m1 / 5
            + 4 * f_p1 / 5
            - 2 * f_p2 / 5
            + 4 * f_p3 / 35
            - f_p4 / 70,
            f_m5 / 252
            - 5 * f_m4 / 126
            + 5 * f_m3 / 28
            - 10 * f_m2 / 21
            + 5 * f_m1 / 6
            + 5 * f_p1 / 6
            - 10 * f_p2 / 21
            + 5 * f_p3 / 28
            - 5 * f_p4 / 126
            + f_p5 / 252,
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"{(half_width * 2)}-point central interpolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    stencil.remove(0)
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil).toSympyExpr()
                    self.assertEqual(expected[half_width], sp.simplify(actual))

        expected = [
            0,
            0,
            2 * f_1 - f_2,
            3 * f_1 - 3 * f_2 + f_3,
            4 * f_1 - 6 * f_2 + 4 * f_3 - f_4,
            5 * f_1 - 10 * f_2 + 10 * f_3 - 5 * f_4 + f_5,
            6 * f_1 - 15 * f_2 + 20 * f_3 - 15 * f_4 + 6 * f_5 - f_6,
            7 * f_1 - 21 * f_2 + 35 * f_3 - 35 * f_4 + 21 * f_5 - 7 * f_6 + f_7,
            8 * f_1
            - 28 * f_2
            + 56 * f_3
            - 70 * f_4
            + 56 * f_5
            - 28 * f_6
            + 8 * f_7
            - f_8,
            9 * f_1
            - 36 * f_2
            + 84 * f_3
            - 126 * f_4
            + 126 * f_5
            - 84 * f_6
            + 36 * f_7
            - 9 * f_8
            + f_9,
        ]
        for shuffle in [False, True]:
            for width in range(2, 10):
                with self.subTest(
                    f"{(width * 2)}-point extrapolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(1, width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil).toSympyExpr()
                    self.assertEqual(expected[width], sp.simplify(actual))

    def test_coefficients(self):
        """
        test suite for interplation.coefficients.
        """

        expected = [
            0,
            [sp.Rational(1, 2), sp.Rational(1, 2)],
            [
                sp.Rational(-1, 6),
                sp.Rational(2, 3),
                sp.Rational(2, 3),
                sp.Rational(-1, 6),
            ],
            [
                sp.Rational(1, 20),
                sp.Rational(-3, 10),
                sp.Rational(3, 4),
                sp.Rational(3, 4),
                sp.Rational(-3, 10),
                sp.Rational(1, 20),
            ],
            [
                sp.Rational(-1, 70),
                sp.Rational(4, 35),
                sp.Rational(-2, 5),
                sp.Rational(4, 5),
                sp.Rational(4, 5),
                sp.Rational(-2, 5),
                sp.Rational(4, 35),
                sp.Rational(-1, 70),
            ],
            [
                sp.Rational(1, 252),
                sp.Rational(-5, 126),
                sp.Rational(5, 28),
                sp.Rational(-10, 21),
                sp.Rational(5, 6),
                sp.Rational(5, 6),
                sp.Rational(-10, 21),
                sp.Rational(5, 28),
                sp.Rational(-5, 126),
                sp.Rational(1, 252),
            ],
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"coefficients of {(half_width * 2)}-point central interpolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    stencil.remove(0)
                    if shuffle:
                        random.shuffle(stencil)
                    actual = coefficients(stencil)
                self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            ([1, 1], 2),
            ([-1, 4, 4, -1], 6),
            ([1, -6, 15, 15, -6, 1], 20),
            ([-1, 8, -28, 56, 56, -28, 8, -1], 70),
            ([1, -10, 45, -120, 210, 210, -120, 45, -10, 1], 252),
        ]
        for half_width in range(1, 6):
            with self.subTest(
                f"coefficients as numerator and denominator of {(half_width * 2)}-point central interpolation"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                stencil.remove(0)
                actual = coefficients(stencil, as_numer_denom=True)
                self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            0,
            [2, -1],
            [3, -3, 1],
            [4, -6, 4, -1],
            [5, -10, 10, -5, 1],
            [6, -15, 20, -15, 6, -1],
            [7, -21, 35, -35, 21, -7, 1],
            [8, -28, 56, -70, 56, -28, 8, -1],
            [9, -36, 84, -126, 126, -84, 36, -9, 1],
        ]
        for shuffle in [False, True]:
            for width in range(2, 10):
                with self.subTest(
                    f"coefficients of {width}-point extrapolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(1, width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = coefficients(stencil)
                    self.assertEqual(expected[width], actual)

        expected = [
            0,
            0,
            ([2, -1], 1),
            ([3, -3, 1], 1),
            ([4, -6, 4, -1], 1),
            ([5, -10, 10, -5, 1], 1),
            ([6, -15, 20, -15, 6, -1], 1),
            ([7, -21, 35, -35, 21, -7, 1], 1),
            ([8, -28, 56, -70, 56, -28, 8, -1], 1),
            ([9, -36, 84, -126, 126, -84, 36, -9, 1], 1),
        ]
        for shuffle in [False, True]:
            for width in range(2, 10):
                with self.subTest(
                    f"coefficients as numerator and denominator of {width}-point extrapolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(1, width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = coefficients(stencil, as_numer_denom=True)
                    self.assertEqual(expected[width], actual)

    def test_truncation_error(self):
        """
        test suite for interplation.truncation_error.
        """
        h = sp.symbols("h")
        f = sp.symbols("f")
        f_1 = sp.symbols("f^(1)")
        f_2 = sp.symbols("f^(2)")
        f_3 = sp.symbols("f^(3)")
        f_4 = sp.symbols("f^(4)")
        f_5 = sp.symbols("f^(5)")
        f_6 = sp.symbols("f^(6)")
        f_7 = sp.symbols("f^(7)")
        f_8 = sp.symbols("f^(8)")
        f_9 = sp.symbols("f^(9)")
        f_10 = sp.symbols("f^(10)")

        expected = [
            0,
            -f_2 * h**2 / 2,
            f_4 * h**4 / 6,
            -f_6 * h**6 / 20,
            f_8 * h**8 / 70,
            -f_10 * h**10 / 252,
        ]
        for half_width in range(1, 6):
            with self.subTest(
                f"truncation error of {(half_width * 2)}-point central interpolation"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                stencil.remove(0)
                actual = truncation_error(stencil)
                self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            0,
            f_2 * h**2,
            -f_3 * h**3,
            f_4 * h**4,
            -f_5 * h**5,
            f_6 * h**6,
            -f_7 * h**7,
            f_8 * h**8,
            -f_9 * h**9,
        ]
        for width in range(2, 10):
            with self.subTest(f"truncation error of {width}-point extrapolation"):
                stencil = [i for i in range(1, width + 1)]
                actual = truncation_error(stencil)
                self.assertEqual(expected[width], actual)


if __name__ == "__main__":
    unittest.main()
