"""Tests for distos.series.taylor_expansion
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import numpy as np
import random

from dictos.series.taylor_expansion import taylor_series, derivative_symbol


class TaylorExpansionTest(unittest.TestCase):
    def test_derivative_symbol(self):
        """
        test suite for taylor_expansion.derivative_symbol.
        """

        for _ in range(10):
            n = random.randint(0, np.iinfo(np.int32).max)
            with self.subTest(n):
                expected = sp.symbols(f"f^({n})")
                actual = derivative_symbol("f", n)
                self.assertEqual(expected, actual)

    def test_taylor_series(self):
        """test suite for taylor_expansion.taylor_series."""

        h = sp.symbols("h")
        a = sp.symbols("a")
        f = sp.symbols("f")
        f_deriv = [derivative_symbol("f", i) for i in range(0, 11)]
        expected = [
            f,
            a * f_deriv[1] * h + f,
            a**2 * f_deriv[2] * h**2 / 2 + a * f_deriv[1] * h + f,
            a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**6 * f_deriv[6] * h**6 / 720
            + a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**7 * f_deriv[7] * h**7 / 5040
            + a**6 * f_deriv[6] * h**6 / 720
            + a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**8 * f_deriv[8] * h**8 / 40320
            + a**7 * f_deriv[7] * h**7 / 5040
            + a**6 * f_deriv[6] * h**6 / 720
            + a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**9 * f_deriv[9] * h**9 / 362880
            + a**8 * f_deriv[8] * h**8 / 40320
            + a**7 * f_deriv[7] * h**7 / 5040
            + a**6 * f_deriv[6] * h**6 / 720
            + a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
            a**10 * f_deriv[10] * h**10 / 3628800
            + a**9 * f_deriv[9] * h**9 / 362880
            + a**8 * f_deriv[8] * h**8 / 40320
            + a**7 * f_deriv[7] * h**7 / 5040
            + a**6 * f_deriv[6] * h**6 / 720
            + a**5 * f_deriv[5] * h**5 / 120
            + a**4 * f_deriv[4] * h**4 / 24
            + a**3 * f_deriv[3] * h**3 / 6
            + a**2 * f_deriv[2] * h**2 / 2
            + a * f_deriv[1] * h
            + f,
        ]

        half_width = 10
        for d in range(-half_width, half_width + 1):
            for i in range(half_width + 1):
                with self.subTest(
                    f"Taylor series around {d}h up to {i} order derivative"
                ):
                    actual = taylor_series(d * h, i)
                    self.assertEqual(expected[i].subs(a, d), actual)


if __name__ == "__main__":
    unittest.main()
