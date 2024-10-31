"""Tests for distos.utilities.utils
"""

import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import numpy as np
import random

from dictos.discrete.stencil import (
    create_coordinate_symbols,
    to_subscript,
)
from dictos.utilities.utils import (
    simplify_coefficients,
    extract_coefficients_as_numer_denom,
    sort_by_subscript,
)
from dictos.linalg.linalg import dot_product
from test.utilities.gen import (
    random_string,
    random_int,
    STENCIL_HALF_WIDTH,
    MAX_SYMBOL_LENGTH,
)
from dictos.discrete.exceptions import TooNarrowError, DuplicatedPointError
from dictos.linalg.exceptions import InconsistentDataSetError


class UtilsTest(unittest.TestCase):
    def test_simplify_coefficients(self):
        """test suite for utils.simplify_coefficients.
        1. it returns integer list when sympy integer list is passed.
        2. it returns integer list when list of sympy symbols is passed.
        3. it returns list of sympy rational when sympy number list is passed.
        4. it returns list of sympy rational when list of sympy symbols is passed.
        5. it returns tuple of integer list and a integer when sympy integer list and as_numer_denom=True are passed.
        6. it returns tuple of integer list and a integer when list of sympy symbols and as_numer_denom=True are passed.
        """

        # subtest 1
        # it returns sympy integers list when sympy integers list is passed.
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                numer = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                expected = [sp.Mul(n, sp.Pow(1, -1)) for n in numer]

                actual = simplify_coefficients(numer)

                self.assertEqual(expected, actual)

        # subtest 2
        # it returns integer list when list of sympy symbols is passed.
        h = sp.symbols(random_string(MAX_SYMBOL_LENGTH))
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                numer = [
                    random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    for _ in range(len_)
                ]
                coef = [n * h for n in numer]
                expected = [n / 1 for n in numer]

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

        # subtest 3
        # it returns list of sympy rational when sympy number list is passed.
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                denom = [sp.Number(random.randint(1, 1000)) for _ in range(len_)]
                numer = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                coef = [sp.Mul(numer[i], sp.Pow(denom[i], -1)) for i in range(len_)]
                expected = [sp.Rational(numer[i], denom[i]) for i in range(len_)]

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

        # subtest 4
        # it returns list of sympy rational when list of sympy symbols is passed.
        h = sp.symbols(random_string(MAX_SYMBOL_LENGTH))
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                denom = [sp.Number(random.randint(1, 1000)) for _ in range(len_)]
                numer = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                coef = [sp.Mul(numer[i], sp.Pow(denom[i], -1)) * h for i in range(len_)]
                expected = [sp.Rational(numer[i], denom[i]) for i in range(len_)]

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

        # subtest 5
        # it returns sympy integers list when sympy integers list is passed.
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                numer = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                expected = ([sp.Mul(n, sp.Pow(1, -1)) for n in numer], 1)

                actual = simplify_coefficients(numer, as_numer_denom=True)

                self.assertEqual(expected, actual)

        # subtest 6
        # it returns tuple of integer list and a integer when list of sympy symbols and as_numer_denom=True are passed.
        h = sp.symbols(random_string(MAX_SYMBOL_LENGTH))
        for len_ in random_int(2, STENCIL_HALF_WIDTH * 2):
            with self.subTest(len_):
                numer = [
                    random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    for _ in range(len_)
                ]

                coef = [n * h for n in numer]
                expected = ([n / 1 for n in numer], 1)

                actual = simplify_coefficients(coef, as_numer_denom=True)

                self.assertEqual(expected, actual)

        # subtest 7
        for coef in [[sp.Number(3 / 4), sp.Number(2 / 3), sp.Number(4 / 5)]]:
            with self.subTest(coef):
                expected = (
                    [45, 40, 48],
                    60,
                )

                actual = simplify_coefficients(coef, as_numer_denom=True)

                self.assertEqual(expected, actual)

    def test_extract_coefficients_as_numer_denom(self):
        """test suite for utils.extract_coefficients_as_numer_denom."""

        x, h, f_0, f_1, f_2, f_3, f_4 = sp.symbols("x h f_0 f_1 f_2 f_3 f_4")
        f_set = [f_0, f_1]
        # 1
        expr = -1.0 * f_0 * (-0.5 * h + x) / h + 1.0 * f_1 * (0.5 * h + x) / h
        with self.subTest(expr):
            expected = ([0.5 * h - 1.0 * x, 0.5 * h + 1.0 * x], [h])
            actual = extract_coefficients_as_numer_denom(expr, f_set)
            self.assertEqual(expected, actual)

        # 2
        expr = 0.5 * f_0 + 0.5 * f_1
        with self.subTest(expr):
            expected = ([0.5, 0.5], [1])
            actual = extract_coefficients_as_numer_denom(expr, f_set)
            self.assertEqual(expected, actual)

        f_set = [f_0, f_1, f_2, f_3, f_4]
        # 3
        expr = (
            f_0 * x * (-2 * h + x) * (-h + x) * (h + x) / (24 * h**4)
            - f_1 * x * (-2 * h + x) * (-h + x) * (2 * h + x) / (6 * h**4)
            + f_2 * (-2 * h + x) * (-h + x) * (h + x) * (2 * h + x) / (4 * h**4)
            - f_3 * x * (-2 * h + x) * (h + x) * (2 * h + x) / (6 * h**4)
            + f_4 * x * (-h + x) * (h + x) * (2 * h + x) / (24 * h**4)
        )
        with self.subTest(expr):
            expected = (
                [
                    2 * h**3 * x - h**2 * x**2 - 2 * h * x**3 + x**4,
                    -16 * h**3 * x + 16 * h**2 * x**2 + 4 * h * x**3 - 4 * x**4,
                    24 * h**4 - 30 * h**2 * x**2 + 6 * x**4,
                    16 * h**3 * x + 16 * h**2 * x**2 - 4 * h * x**3 - 4 * x**4,
                    -2 * h**3 * x - h**2 * x**2 + 2 * h * x**3 + x**4,
                ],
                [24 * h**4],
            )
            actual = extract_coefficients_as_numer_denom(expr, f_set)
            self.assertEqual(expected, actual)

        # 4
        expr = -f_0 / 6 + 2 * f_1 / 3 + 2 * f_2 / 3 - f_3 / 6
        with self.subTest(expr):
            expected = ([-1, 4, 4, -1], [6])
            actual = extract_coefficients_as_numer_denom(expr, f_set)
            self.assertEqual(expected, actual)

    def test_utils_exception(self):
        """test suite for exception in utils"""

        stencil = [0]
        with self.subTest("create_coordinate_symbols with too narrow stencil"):
            with self.assertRaises(TooNarrowError):
                create_coordinate_symbols(stencil)

        stencil = [1, 1, 2, 3, 4]
        with self.subTest("create_coordinate_symbols with invalid stencil"):
            with self.assertRaises(DuplicatedPointError):
                create_coordinate_symbols(stencil)

    def test_utils_sort_by_subscript(self):
        """test suite for utils.sort_by_subscript."""

        for half_width in range(1, 11):
            with self.subTest(f"get subscript {(half_width * 2 + 1)}-point stencil"):
                stencil = [to_subscript(i) for i in range(-half_width, half_width + 1)]
                sym_str = "".join(["f" + "_{" + s + "}" + " " for s in stencil])
                f_set = sp.symbols(sym_str)
                expected = dot_product(
                    [2 for _ in range(len(f_set))], f_set, evaluate=False
                )

                num = random_int(-half_width, half_width)
                stencil = [to_subscript(i) for i in num]
                sym_str = "".join(["f" + "_{" + s + "}" + " " for s in stencil])
                f_set = sp.symbols(sym_str)
                eq = dot_product(f_set, [2 for _ in range(len(f_set))], evaluate=False)
                actual = sort_by_subscript(eq)

                ac_str = str(actual)
                ex_str = str(expected)
                self.assertEqual(ex_str, ac_str)

                stencil = [to_subscript(i) for i in range(-half_width, half_width + 1)]
                sym_str = "".join(["f" + "_{" + s + "}" + " " for s in stencil])
                f_set = sp.symbols(sym_str)
                expected = sp.Add(*f_set, evaluate=False)

                num = random_int(-half_width, half_width)
                stencil = [to_subscript(i) for i in num]
                sym_str = "".join(["f" + "_{" + s + "}" + " " for s in stencil])
                f_set = sp.symbols(sym_str)
                eq = sp.Add(*f_set, evaluate=False)
                actual = sort_by_subscript(eq)

                ac_str = str(actual)
                ex_str = str(expected)
                self.assertEqual(ex_str, ac_str)


if __name__ == "__main__":
    unittest.main()
