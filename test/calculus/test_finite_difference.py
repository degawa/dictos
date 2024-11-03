"""Tests for distos.calculus.finite_difference
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.calculus.finite_difference import (
    equation,
    coefficients,
    truncation_error,
    generate,
)


class FiniteDifferenceTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for finite_difference.equation.
        """

        h = sp.Symbol("h")
        f_0 = sp.Symbol("f_{0}")
        f_1 = sp.Symbol("f_{1}")
        f_2 = sp.Symbol("f_{2}")
        f_3 = sp.Symbol("f_{3}")
        f_4 = sp.Symbol("f_{4}")
        f_5 = sp.Symbol("f_{5}")
        f_6 = sp.Symbol("f_{6}")
        f_7 = sp.Symbol("f_{7}")
        f_8 = sp.Symbol("f_{8}")
        f_9 = sp.Symbol("f_{9}")
        f_10 = sp.Symbol("f_{10}")

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
            (-f_m1 + f_p1) / (2 * h),  # type: ignore
            (f_m2 - 8 * f_m1 + 8 * f_p1 - f_p2) / (12 * h),  # type: ignore
            (-f_m3 + 9 * f_m2 - 45 * f_m1 + 45 * f_p1 - 9 * f_p2 + f_p3) / (60 * h),  # type: ignore
            (
                3 * f_m4  # type: ignore
                - 32 * f_m3  # type: ignore
                + 168 * f_m2  # type: ignore
                - 672 * f_m1  # type: ignore
                + 672 * f_p1  # type: ignore
                - 168 * f_p2  # type: ignore
                + 32 * f_p3  # type: ignore
                - 3 * f_p4  # type: ignore
            )
            / (840 * h),  # type: ignore
            (
                -2 * f_m5  # type: ignore
                + 25 * f_m4  # type: ignore
                - 150 * f_m3  # type: ignore
                + 600 * f_m2  # type: ignore
                - 2100 * f_m1  # type: ignore
                + 2100 * f_p1  # type: ignore
                - 600 * f_p2  # type: ignore
                + 150 * f_p3  # type: ignore
                - 25 * f_p4  # type: ignore
                + 2 * f_p5  # type: ignore
            )
            / (2520 * h),  # type: ignore
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"{(half_width * 2 + 1)}-point central difference for 1st derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 1).toSympyExpr()
                    self.assertEqual(expected[half_width], sp.simplify(actual))

        expected = [
            0,
            (f_m1 - 2 * f_0 + f_p1) / h**2,  # type: ignore
            (-f_m2 + 16 * f_m1 - 30 * f_0 + 16 * f_p1 - f_p2) / (12 * h**2),  # type: ignore
            (
                2 * f_m3  # type: ignore
                - 27 * f_m2  # type: ignore
                + 270 * f_m1  # type: ignore
                - 490 * f_0  # type: ignore
                + 270 * f_p1  # type: ignore
                - 27 * f_p2  # type: ignore
                + 2 * f_p3  # type: ignore
            )
            / (180 * h**2),  # type: ignore
            (
                -9 * f_m4  # type: ignore
                + 128 * f_m3  # type: ignore
                - 1008 * f_m2  # type: ignore
                + 8064 * f_m1  # type: ignore
                - 14350 * f_0  # type: ignore
                + 8064 * f_p1  # type: ignore
                - 1008 * f_p2  # type: ignore
                + 128 * f_p3  # type: ignore
                - 9 * f_p4  # type: ignore
            )
            / (5040 * h**2),  # type: ignore
            (
                8 * f_m5  # type: ignore
                - 125 * f_m4  # type: ignore
                + 1000 * f_m3  # type: ignore
                - 6000 * f_m2  # type: ignore
                + 42000 * f_m1  # type: ignore
                - 73766 * f_0  # type: ignore
                + 42000 * f_p1  # type: ignore
                - 6000 * f_p2  # type: ignore
                + 1000 * f_p3  # type: ignore
                - 125 * f_p4  # type: ignore
                + 8 * f_p5  # type: ignore
            )
            / (25200 * h**2),  # type: ignore
        ]
        for shuffle in [False, True]:
            for half_width in range(1, 6):
                with self.subTest(
                    f"{(half_width * 2 + 1)}-point central difference for 2nd derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(-half_width, half_width + 1)]
                    if shuffle:
                        random.shuffle(stencil)
                    actual = equation(stencil, 2).toSympyExpr()
                    self.assertEqual(expected[half_width], sp.simplify(actual))

        expected = [
            0,
            (-f_0 + f_1) / h,  # type: ignore
            (-3 * f_0 + 4 * f_1 - f_2) / (2 * h),  # type: ignore
            (-11 * f_0 + 18 * f_1 - 9 * f_2 + 2 * f_3) / (6 * h),  # type: ignore
            (-25 * f_0 + 48 * f_1 - 36 * f_2 + 16 * f_3 - 3 * f_4) / (12 * h),  # type: ignore
            (-137 * f_0 + 300 * f_1 - 300 * f_2 + 200 * f_3 - 75 * f_4 + 12 * f_5)  # type: ignore
            / (60 * h),  # type: ignore
            (
                -147 * f_0  # type: ignore
                + 360 * f_1  # type: ignore
                - 450 * f_2  # type: ignore
                + 400 * f_3  # type: ignore
                - 225 * f_4  # type: ignore
                + 72 * f_5  # type: ignore
                - 10 * f_6  # type: ignore
            )
            / (60 * h),  # type: ignore
            (
                -1089 * f_0  # type: ignore
                + 2940 * f_1  # type: ignore
                - 4410 * f_2  # type: ignore
                + 4900 * f_3  # type: ignore
                - 3675 * f_4  # type: ignore
                + 1764 * f_5  # type: ignore
                - 490 * f_6  # type: ignore
                + 60 * f_7  # type: ignore
            )
            / (420 * h),  # type: ignore
            (
                -2283 * f_0  # type: ignore
                + 6720 * f_1  # type: ignore
                - 11760 * f_2  # type: ignore
                + 15680 * f_3  # type: ignore
                - 14700 * f_4  # type: ignore
                + 9408 * f_5  # type: ignore
                - 3920 * f_6  # type: ignore
                + 960 * f_7  # type: ignore
                - 105 * f_8  # type: ignore
            )
            / (840 * h),  # type: ignore
            (
                -7129 * f_0  # type: ignore
                + 22680 * f_1  # type: ignore
                - 45360 * f_2  # type: ignore
                + 70560 * f_3  # type: ignore
                - 79380 * f_4  # type: ignore
                + 63504 * f_5  # type: ignore
                - 35280 * f_6  # type: ignore
                + 12960 * f_7  # type: ignore
                - 2835 * f_8  # type: ignore
                + 280 * f_9  # type: ignore
            )
            / (2520 * h),  # type: ignore
            (
                -7381 * f_0  # type: ignore
                + 25200 * f_1  # type: ignore
                - 56700 * f_2  # type: ignore
                + 100800 * f_3  # type: ignore
                - 132300 * f_4  # type: ignore
                + 127008 * f_5  # type: ignore
                - 88200 * f_6  # type: ignore
                + 43200 * f_7  # type: ignore
                - 14175 * f_8  # type: ignore
                + 2800 * f_9  # type: ignore
                - 252 * f_10  # type: ignore
            )
            / (2520 * h),  # type: ignore
        ]
        for shuffle in [False, True]:
            for width in range(1, 11):
                with self.subTest(
                    f"{(width + 1)}-point forward difference for 1st derivative, stencil shuffle = {shuffle}"
                ):
                    stencil = [i for i in range(width + 1)]
                    random.shuffle(stencil)
                    actual = equation(stencil, 1).toSympyExpr()
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

    def test_generate(self):
        """
        test suite for finite_difference.generate.
        """

        grid_type = "regular"

        acc = 2
        deriv_stencils = [
            (1, [-1, 0, 1]),
            (2, [-1, 0, 1]),
            (3, [-2, -1, 0, 1, 2]),
            (4, [-2, -1, 0, 1, 2]),
            (5, [-3, -2, -1, 0, 1, 2, 3]),
            (6, [-3, -2, -1, 0, 1, 2, 3]),
            (7, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (8, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (9, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (10, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
        ]

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference equation for {deriv}-derivative on regular grid"
            ):

                expected = equation(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=True
                )
                self.assertEqual(expected, actual)

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference coefficients for {deriv}-derivative on regular grid"
            ):

                expected = coefficients(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=False
                )
                self.assertEqual(expected, actual)

        acc = 4
        deriv_stencils = [
            (1, [-2, -1, 0, 1, 2]),
            (2, [-2, -1, 0, 1, 2]),
            (3, [-3, -2, -1, 0, 1, 2, 3]),
            (4, [-3, -2, -1, 0, 1, 2, 3]),
            (5, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (6, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (7, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (8, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (9, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (10, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
        ]

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference equation for {deriv}-derivative on regular grid"
            ):

                expected = equation(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=True
                )
                self.assertEqual(expected, actual)

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference coefficients for {deriv}-derivative on regular grid"
            ):

                expected = coefficients(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=False
                )
                self.assertEqual(expected, actual)

        acc = 6
        deriv_stencils = [
            (1, [-3, -2, -1, 0, 1, 2, 3]),
            (2, [-3, -2, -1, 0, 1, 2, 3]),
            (3, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (4, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (5, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (6, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (7, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (8, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (9, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
            (10, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
        ]

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference equation for {deriv}-derivative on regular grid"
            ):

                expected = equation(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=True
                )
                self.assertEqual(expected, actual)

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference coefficients for {deriv}-derivative on regular grid"
            ):

                expected = coefficients(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=False
                )
                self.assertEqual(expected, actual)

        acc = 8
        deriv_stencils = [
            (1, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (2, [-4, -3, -2, -1, 0, 1, 2, 3, 4]),
            (3, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (4, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (5, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (6, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (7, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
            (8, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
            (9, [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (10, [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
        ]

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference equation for {deriv}-derivative on regular grid"
            ):

                expected = equation(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=True
                )
                self.assertEqual(expected, actual)

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference coefficients for {deriv}-derivative on regular grid"
            ):

                expected = coefficients(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=False
                )
                self.assertEqual(expected, actual)

        acc = 10
        deriv_stencils = [
            (1, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (2, [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]),
            (3, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (4, [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]),
            (5, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
            (6, [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]),
            (7, [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (8, [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (9, [-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            (10, [-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ]

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference equation for {deriv}-derivative on regular grid"
            ):

                expected = equation(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=True
                )
                self.assertEqual(expected, actual)

        for deriv, stencil in deriv_stencils:
            with self.subTest(
                f"generate {acc}-order central difference coefficients for {deriv}-derivative on regular grid"
            ):

                expected = coefficients(stencil, deriv=deriv)
                actual = generate(
                    deriv=deriv, acc=acc, grid_type=grid_type, as_equation=False
                )
                self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
