"""Tests for distos.finite_difference
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.calculus.finite_difference import equation, coefficients, truncation_error


class FiniteDifferenceTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for finite_difference.equation.
        """

        h = sp.symbols("h")
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
            (-f_m1 + f_p1) / (2 * h),
            (f_m2 - 8 * f_m1 + 8 * f_p1 - f_p2) / (12 * h),
            (-f_m3 + 9 * f_m2 - 45 * f_m1 + 45 * f_p1 - 9 * f_p2 + f_p3) / (60 * h),
            (
                3 * f_m4
                - 32 * f_m3
                + 168 * f_m2
                - 672 * f_m1
                + 672 * f_p1
                - 168 * f_p2
                + 32 * f_p3
                - 3 * f_p4
            )
            / (840 * h),
            (
                -2 * f_m5
                + 25 * f_m4
                - 150 * f_m3
                + 600 * f_m2
                - 2100 * f_m1
                + 2100 * f_p1
                - 600 * f_p2
                + 150 * f_p3
                - 25 * f_p4
                + 2 * f_p5
            )
            / (2520 * h),
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"{(half_width * 2 + 1)}-point central difference for 1st derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 1)
                    self.assertEqual(expected[half_width], sp.simplify(actual))

        expected = [
            0,
            (f_m1 - 2 * f_0 + f_p1) / h**2,
            (-f_m2 + 16 * f_m1 - 30 * f_0 + 16 * f_p1 - f_p2) / (12 * h**2),
            (
                2 * f_m3
                - 27 * f_m2
                + 270 * f_m1
                - 490 * f_0
                + 270 * f_p1
                - 27 * f_p2
                + 2 * f_p3
            )
            / (180 * h**2),
            (
                -9 * f_m4
                + 128 * f_m3
                - 1008 * f_m2
                + 8064 * f_m1
                - 14350 * f_0
                + 8064 * f_p1
                - 1008 * f_p2
                + 128 * f_p3
                - 9 * f_p4
            )
            / (5040 * h**2),
            (
                8 * f_m5
                - 125 * f_m4
                + 1000 * f_m3
                - 6000 * f_m2
                + 42000 * f_m1
                - 73766 * f_0
                + 42000 * f_p1
                - 6000 * f_p2
                + 1000 * f_p3
                - 125 * f_p4
                + 8 * f_p5
            )
            / (25200 * h**2),
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"{(half_width * 2 + 1)}-point central difference for 2nd derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 2)
                    self.assertEqual(expected[half_width], sp.simplify(actual))

        expected = [
            0,
            (-f_0 + f_1) / h,
            (-3 * f_0 + 4 * f_1 - f_2) / (2 * h),
            (-11 * f_0 + 18 * f_1 - 9 * f_2 + 2 * f_3) / (6 * h),
            (-25 * f_0 + 48 * f_1 - 36 * f_2 + 16 * f_3 - 3 * f_4) / (12 * h),
            (-137 * f_0 + 300 * f_1 - 300 * f_2 + 200 * f_3 - 75 * f_4 + 12 * f_5)
            / (60 * h),
            (
                -147 * f_0
                + 360 * f_1
                - 450 * f_2
                + 400 * f_3
                - 225 * f_4
                + 72 * f_5
                - 10 * f_6
            )
            / (60 * h),
            (
                -1089 * f_0
                + 2940 * f_1
                - 4410 * f_2
                + 4900 * f_3
                - 3675 * f_4
                + 1764 * f_5
                - 490 * f_6
                + 60 * f_7
            )
            / (420 * h),
            (
                -2283 * f_0
                + 6720 * f_1
                - 11760 * f_2
                + 15680 * f_3
                - 14700 * f_4
                + 9408 * f_5
                - 3920 * f_6
                + 960 * f_7
                - 105 * f_8
            )
            / (840 * h),
            (
                -7129 * f_0
                + 22680 * f_1
                - 45360 * f_2
                + 70560 * f_3
                - 79380 * f_4
                + 63504 * f_5
                - 35280 * f_6
                + 12960 * f_7
                - 2835 * f_8
                + 280 * f_9
            )
            / (2520 * h),
            (
                -7381 * f_0
                + 25200 * f_1
                - 56700 * f_2
                + 100800 * f_3
                - 132300 * f_4
                + 127008 * f_5
                - 88200 * f_6
                + 43200 * f_7
                - 14175 * f_8
                + 2800 * f_9
                - 252 * f_10
            )
            / (2520 * h),
        ]
        for shuffle in [False, True]:
            for width in range(1, 11):
                with self.subTest(
                    f"{(width + 1)}-point forward difference for 1st derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(width + 1)]
                    random.shuffle(stencil)
                    actual = equation(stencil, 1)
                    self.assertEqual(expected[width], sp.simplify(actual))

    def test_coefficients(self):
        """
        test suite for finite_difference.coefficients.
        """

        expected = [
            0,
            [sp.Rational(-1, 2), 0, sp.Rational(1, 2)],
            [
                sp.Rational(1, 12),
                sp.Rational(-2, 3),
                0,
                sp.Rational(2, 3),
                sp.Rational(-1, 12),
            ],
            [
                sp.Rational(-1, 60),
                sp.Rational(3, 20),
                sp.Rational(-3, 4),
                0,
                sp.Rational(3, 4),
                sp.Rational(-3, 20),
                sp.Rational(1, 60),
            ],
            [
                sp.Rational(1, 280),
                sp.Rational(-4, 105),
                sp.Rational(1, 5),
                sp.Rational(-4, 5),
                0,
                sp.Rational(4, 5),
                sp.Rational(-1, 5),
                sp.Rational(4, 105),
                sp.Rational(-1, 280),
            ],
            [
                sp.Rational(-1, 1260),
                sp.Rational(5, 504),
                sp.Rational(-5, 84),
                sp.Rational(5, 21),
                sp.Rational(-5, 6),
                0,
                sp.Rational(5, 6),
                sp.Rational(-5, 21),
                sp.Rational(5, 84),
                sp.Rational(-5, 504),
                sp.Rational(1, 1260),
            ],
        ]

        for half_width in range(1, 6):
            with self.subTest(
                f"coefficents of {(half_width * 2 + 1)}-point central difference for 1st derivative"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = coefficients(stencil, 1)
                self.assertEqual(expected[half_width], actual)

        for half_width in range(1, 6):
            with self.subTest(
                f"coefficents of {(half_width * 2 + 1)}-point central difference for 1st derivative"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                random.shuffle(stencil)
                actual = coefficients(stencil, 1)
                self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            ([-1, 0, 1], 2),
            ([1, -8, 0, 8, -1], 12),
            ([-1, 9, -45, 0, 45, -9, 1], 60),
            ([3, -32, 168, -672, 0, 672, -168, 32, -3], 840),
            ([-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2], 2520),
        ]
        for half_width in range(1, 6):
            with self.subTest(
                f"numerator and denominator of coefficients of {(half_width * 2 + 1)}-point central difference for 1st derivative"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = coefficients(stencil, 1, as_numer_denom=True)
                self.assertEqual(expected[half_width], actual)

        for half_width in range(1, 6):
            with self.subTest(
                f"numerator and denominator of coefficients of {(half_width * 2 + 1)}-point central difference for 1st derivative"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                random.shuffle(stencil)
                actual = coefficients(stencil, 1, as_numer_denom=True)
                self.assertEqual(expected[half_width], actual)

    def test_truncation_error(self):
        """
        test suite for finite_difference.truncation_error.
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
        f_11 = sp.symbols("f^(11)")

        expected = [
            0,
            -f_3 * h**2 / 6,
            f_5 * h**4 / 30,
            -f_7 * h**6 / 140,
            f_9 * h**8 / 630,
            -f_11 * h**10 / 2772,
        ]
        for half_width in range(1, 6):
            with self.subTest(
                f"truncation error of {(half_width * 2 + 1)}-point central difference for 1st derivative"
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = truncation_error(stencil, 1)
                self.assertEqual(expected[half_width], actual)


if __name__ == "__main__":
    unittest.main()
