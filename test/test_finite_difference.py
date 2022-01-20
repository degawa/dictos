"""Tests for distos.finite_difference
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.finite_difference import equation, coefficients, truncation_error


class FiniteDifferenceTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for finite_difference.equation.
        """

        h = sp.symbols("h")
        f_0 = sp.symbols("f_0")
        f_1 = sp.symbols("f_1")
        f_2 = sp.symbols("f_2")
        f_3 = sp.symbols("f_3")
        f_4 = sp.symbols("f_4")
        f_5 = sp.symbols("f_5")
        f_6 = sp.symbols("f_6")
        f_7 = sp.symbols("f_7")
        f_8 = sp.symbols("f_8")
        f_9 = sp.symbols("f_9")
        f_10 = sp.symbols("f_10")

        expected = [
            0,
            (-f_0 + f_2) / (2 * h),
            (f_0 - 8 * f_1 + 8 * f_3 - f_4) / (12 * h),
            (-f_0 + 9 * f_1 - 45 * f_2 + 45 * f_4 - 9 * f_5 + f_6) / (60 * h),
            (
                3 * f_0
                - 32 * f_1
                + 168 * f_2
                - 672 * f_3
                + 672 * f_5
                - 168 * f_6
                + 32 * f_7
                - 3 * f_8
            )
            / (840 * h),
            (
                -2 * f_0
                + 25 * f_1
                + 2 * f_10
                - 150 * f_2
                + 600 * f_3
                - 2100 * f_4
                + 2100 * f_6
                - 600 * f_7
                + 150 * f_8
                - 25 * f_9
            )
            / (2520 * h),
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    "%d-point central difference for 1st derivative, stencil shuffle = %r"
                    % (half_width * 2 + 1, shuffle)
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 1)
                    self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            (f_0 - 2 * f_1 + f_2) / h ** 2,
            (-f_0 + 16 * f_1 - 30 * f_2 + 16 * f_3 - f_4) / (12 * h ** 2),
            (
                2 * f_0
                - 27 * f_1
                + 270 * f_2
                - 490 * f_3
                + 270 * f_4
                - 27 * f_5
                + 2 * f_6
            )
            / (180 * h ** 2),
            (
                -9 * f_0
                + 128 * f_1
                - 1008 * f_2
                + 8064 * f_3
                - 14350 * f_4
                + 8064 * f_5
                - 1008 * f_6
                + 128 * f_7
                - 9 * f_8
            )
            / (5040 * h ** 2),
            (
                8 * f_0
                - 125 * f_1
                + 8 * f_10
                + 1000 * f_2
                - 6000 * f_3
                + 42000 * f_4
                - 73766 * f_5
                + 42000 * f_6
                - 6000 * f_7
                + 1000 * f_8
                - 125 * f_9
            )
            / (25200 * h ** 2),
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    "%d-point central difference for 2nd derivative, stencil shuffle = %r"
                    % (half_width * 2 + 1, shuffle),
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 2)
                    self.assertEqual(expected[half_width], actual)

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
                - 252 * f_10
                - 56700 * f_2
                + 100800 * f_3
                - 132300 * f_4
                + 127008 * f_5
                - 88200 * f_6
                + 43200 * f_7
                - 14175 * f_8
                + 2800 * f_9
            )
            / (2520 * h),
        ]
        for shuffle in [False, True]:
            for width in range(1, 11):
                with self.subTest(
                    "%d-point forward difference for 1st derivative, stencil shuffle = %r"
                    % (width + 1, shuffle)
                ):
                    stencil = [i for i in range(width + 1)]
                    random.shuffle(stencil)
                    actual = equation(stencil, 1)
                    self.assertEqual(expected[width], actual)

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
                "coefficents of %d-point central difference for 1st derivative"
                % (half_width * 2 + 1)
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = coefficients(stencil, 1)
                self.assertEqual(expected[half_width], actual)

        for half_width in range(1, 6):
            with self.subTest(
                "coefficents of %d-point central difference for 1st derivative"
                % (half_width * 2 + 1)
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
                "numerator and denominator of coefficients of %d-point central difference for 1st derivative"
                % (half_width * 2 + 1)
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = coefficients(stencil, 1, as_numer_denom=True)
                self.assertEqual(expected[half_width], actual)

        for half_width in range(1, 6):
            with self.subTest(
                "numerator and denominator of coefficients of %d-point central difference for 1st derivative"
                % (half_width * 2 + 1)
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
            -f_3 * h ** 2 / 6,
            f_5 * h ** 4 / 30,
            -f_7 * h ** 6 / 140,
            f_9 * h ** 8 / 630,
            -f_11 * h ** 10 / 2772,
        ]
        for half_width in range(1, 6):
            with self.subTest(
                "truncation error of %d-point central difference for 1st derivative"
                % (half_width * 2 + 1)
            ):
                stencil = [i for i in range(-half_width, half_width + 1)]
                actual = truncation_error(stencil, 1)
                self.assertEqual(expected[half_width], actual)


if __name__ == "__main__":
    unittest.main()
