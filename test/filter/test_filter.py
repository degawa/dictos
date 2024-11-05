"""Tests for distos.filter.filter
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.filter.filter import generate
from dictos.linalg.linalg import scale


class FilterTest(unittest.TestCase):
    def test_generate(self):
        """
        test suite for filter.generate.
        """

        EXPECTED_FILTER_COEFFICIENTS = {
            2: ([1, 2, 1], 4),
            4: ([-1, 4, 10, 4, -1], 16),
            6: ([1, -6, 15, 44, 15, -6, 1], 64),
            8: ([-1, 8, -28, 56, 186, 56, -28, 8, -1], 256),
            10: ([1, -10, 45, -120, 210, 772, 210, -120, 45, -10, 1], 1024),
        }

        for acc, expected in EXPECTED_FILTER_COEFFICIENTS.items():
            with self.subTest(
                f"generate {acc}-order central filter coefficients as_numer_denom on regular grid"
            ):

                actual = generate(acc=acc, as_numer_denom=True)
                self.assertEqual(expected, actual)

        for acc, expected_nmr_dnm in EXPECTED_FILTER_COEFFICIENTS.items():
            with self.subTest(
                f"generate {acc}-order central filter coefficients on regular grid"
            ):

                numer, denom = expected_nmr_dnm
                expected = scale(numer, sp.Pow(denom, -1))
                actual = generate(acc=acc)
                self.assertEqual(expected, actual)

        f_0 = sp.Symbol("f_{0}")
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
        EXPECTED_FILTER_EQUATION = {
            2: sp.nsimplify((1 * f_m1 + 2 * f_0 + 1 * f_p1) * sp.Pow(4, -1)),  # type: ignore
            4: sp.nsimplify((-1 * f_m2 + 4 * f_m1 + 10 * f_0 + 4 * f_p1 - 1 * f_p2) * sp.Pow(16, -1)),  # type: ignore
            6: sp.nsimplify(
                (
                    1 * f_m3  # type: ignore
                    - 6 * f_m2  # type: ignore
                    + 15 * f_m1  # type: ignore
                    + 44 * f_0  # type: ignore
                    + 15 * f_p1  # type: ignore
                    + -6 * f_p2  # type: ignore
                    + 1 * f_p3  # type: ignore
                )
                * sp.Pow(64, -1)
            ),
            8: sp.nsimplify(
                (
                    -1 * f_m4  # type: ignore
                    + 8 * f_m3  # type: ignore
                    - 28 * f_m2  # type: ignore
                    + 56 * f_m1  # type: ignore
                    + 186 * f_0  # type: ignore
                    + 56 * f_p1  # type: ignore
                    + -28 * f_p2  # type: ignore
                    + 8 * f_p3  # type: ignore
                    - 1 * f_p4  # type: ignore
                )
                * sp.Pow(256, -1)
            ),
            10: sp.nsimplify(
                (
                    1 * f_m5  # type: ignore
                    - 10 * f_m4  # type: ignore
                    + 45 * f_m3  # type: ignore
                    - 120 * f_m2  # type: ignore
                    + 210 * f_m1  # type: ignore
                    + 772 * f_0  # type: ignore
                    + 210 * f_p1  # type: ignore
                    - 120 * f_p2  # type: ignore
                    + 45 * f_p3  # type: ignore
                    - 10 * f_p4  # type: ignore
                    + 1 * f_p5  # type: ignore
                )
                * sp.Pow(1024, -1)
            ),
        }
        for acc, expected in EXPECTED_FILTER_EQUATION.items():
            with self.subTest(
                f"generate {acc}-order central filter equation on regular grid"
            ):

                actual = sp.simplify(generate(acc=acc, as_equation=True).toSympyExpr())
                # could not compare directly Expr and sympy.Expr
                # so compare between sympy.Expr and converted and simplified form
                self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
