"""Tests for distos.interplation
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.interpolation import equation, coefficients, truncation_error


class InterpolationTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for interplation.equation.
        """

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
            f_0 / 2 + f_1 / 2,
            -f_0 / 6 + 2 * f_1 / 3 + 2 * f_2 / 3 - f_3 / 6,
            f_0 / 20
            - 3 * f_1 / 10
            + 3 * f_2 / 4
            + 3 * f_3 / 4
            - 3 * f_4 / 10
            + f_5 / 20,
            -f_0 / 70
            + 4 * f_1 / 35
            - 2 * f_2 / 5
            + 4 * f_3 / 5
            + 4 * f_4 / 5
            - 2 * f_5 / 5
            + 4 * f_6 / 35
            - f_7 / 70,
            f_0 / 252
            - 5 * f_1 / 126
            + 5 * f_2 / 28
            - 10 * f_3 / 21
            + 5 * f_4 / 6
            + 5 * f_5 / 6
            - 10 * f_6 / 21
            + 5 * f_7 / 28
            - 5 * f_8 / 126
            + f_9 / 252,
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
                    actual = equation(stencil)
                    self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            0,
            2 * f_0 - f_1,
            3 * f_0 - 3 * f_1 + f_2,
            4 * f_0 - 6 * f_1 + 4 * f_2 - f_3,
            5 * f_0 - 10 * f_1 + 10 * f_2 - 5 * f_3 + f_4,
            6 * f_0 - 15 * f_1 + 20 * f_2 - 15 * f_3 + 6 * f_4 - f_5,
            7 * f_0 - 21 * f_1 + 35 * f_2 - 35 * f_3 + 21 * f_4 - 7 * f_5 + f_6,
            8 * f_0
            - 28 * f_1
            + 56 * f_2
            - 70 * f_3
            + 56 * f_4
            - 28 * f_5
            + 8 * f_6
            - f_7,
            9 * f_0
            - 36 * f_1
            + 84 * f_2
            - 126 * f_3
            + 126 * f_4
            - 84 * f_5
            + 36 * f_6
            - 9 * f_7
            + f_8,
        ]
        for shuffle in [False, True]:
            for width in range(2, 10):
                with self.subTest(
                    f"{(width * 2)}-point extrapolation, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(1, width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil)
                    self.assertEqual(expected[width], actual)

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
            -f_2 * h ** 2 / 2,
            f_4 * h ** 4 / 6,
            -f_6 * h ** 6 / 20,
            f_8 * h ** 8 / 70,
            -f_10 * h ** 10 / 252,
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
            f_2 * h ** 2,
            -f_3 * h ** 3,
            f_4 * h ** 4,
            -f_5 * h ** 5,
            f_6 * h ** 6,
            -f_7 * h ** 7,
            f_8 * h ** 8,
            -f_9 * h ** 9,
        ]
        for width in range(2, 10):
            with self.subTest(f"truncation error of {width}-point extrapolation"):
                stencil = [i for i in range(1, width + 1)]
                actual = truncation_error(stencil)
                self.assertEqual(expected[width], actual)


if __name__ == "__main__":
    unittest.main()
