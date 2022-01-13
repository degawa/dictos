"""Tests for distos.utils
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import numpy as np
import random
import string

from dictos.utils import (
    DEFAULT_INTERVAL,
    DEFAULT_FUNCTION,
    create_coordinate_symbols,
    create_function_symbols,
    simplify_coefficients,
    dot_product,
    div,
)

_stencil_half_width = 20  # up to 20th order accuracy
_max_symbol_length = 10


def random_string(len):
    """generate n-length random string

    Args:
        len (integer): length of string
    """

    return "".join(random.choices(string.ascii_letters + string.digits, k=len))


def random_int(min, max, exclude=None):
    """generate random list of integers from min to max.
    Integers listed in exclude are excluded

    Args:
        min (integer): minimum value of integer in the list.
        max (integer): maximum value of integer in the list.
        exclude (list of integer): integers to be excluded from the list.
            Defaults to None.
    """
    num = list(range(min, max + 1))

    if exclude is not None:
        for i in exclude:
            num.remove(i)

    random.shuffle(num)

    return num


class UtilsTest(unittest.TestCase):
    def test_create_coordinate_symbols(self):
        """
        test suite for utils.create_coordinate_symbols.
        1. it returns [n*h] when [n] is passed.
        2. it returns [a*h b*h c*h ...] when [a b c ...] is passed.
        3. it returns [a*dx b*dx c*dx ...] when [a b c ...] and 'dx' are passed.
        4. it raise error when empty list is passed.
        5. it raise error when at least a number in the stencil appears more than once.
        """

        # subtest 1
        # it returns [n*h] when [n] is passed.
        num = random_int(-_stencil_half_width, _stencil_half_width, exclude=[0])
        for n in num:
            with self.subTest(n):
                stencil = [n]
                expected = create_coordinate_symbols(stencil)

                h = sp.symbols(DEFAULT_INTERVAL)
                acctual = [n * h]
                self.assertEqual(expected, acctual)

        # subtest 2
        # it returns [a*h b*h c*h ...] when [a b c ...] is passed.
        num = random_int(1, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil)

                h = sp.symbols(DEFAULT_INTERVAL)
                acctual = [i * h for i in range(n)]
                self.assertEqual(expected, acctual)

        # uncomment test after implement error raising when len(stencil)==0
        # with self.subTest():
        #     with self.assertRaises(Exception):
        #         create_coordinate_symbols([])

        # uncomment test after implement error raising when at least a number in the stencil appears more than once.
        # with self.subTest():
        #     with self.assertRaises(Exception):
        #         create_coordinate_symbols([1, 1])

        # subtest 3
        # it returns [a*dx b*dx c*dx ...] when [a b c ...] and 'dx' are passed.
        num = random_int(1, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                interval = random_string(random.randint(1, _max_symbol_length))
                h = sp.symbols(interval)

                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil, interval=interval)

                acctual = [i * h for i in range(n)]
                self.assertEqual(expected, acctual)

    def test_create_function_symbols(self):
        """
        test suite for utils.create_function_symbols.
        1. it returns [f_{0}] when [n*h] is passed.
        2. it returns [f_{0}, f_{1}, f_{2}, ...] when [a*h, b*h, c*h, ...] is passed.
        3. it returns [g_{0}] when [n*h] and 'g' are passed.
        4. it returns [g_{0}, g_{1}, g_{2}, ...] when [a*h, b*h, c*h, ...] and 'g' are passed.
        5. it returns [f_{n}] when [n*h] and same_subscripts_as_stencil=True are passed.
        6. it returns [f_{a}, f_{b}, f_{c}, ...] when [a*h, b*h, c*h, ...] and same_subscripts_as_stencil=True are passed.
        7. it returns [g_{n}] when [n*h], 'g', and same_subscripts_as_stencil=True are passed.
        8. it returns [g_{a}, g_{b}, g_{c}, ...] when [a*h, b*h, c*h, ...], 'g', and same_subscripts_as_stencil=True are passed.
        9. it raise error when empty list is passed.
        """

        # subtest 1
        # it returns [f_{0}] when [n*h] is passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                x = create_coordinate_symbols([n])
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                acctual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, acctual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2])
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                acctual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, acctual)

        # subtest 2
        # it returns [f_{0}, f_{1}, f_{2}, ...] when [a*h, b*h, c*h, ...] is passed.
        num = random_int(1, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                acctual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, acctual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                acctual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, acctual)

        # subtest 3
        # it returns [g_{0}] when [n*h] and 'g' are passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                x = create_coordinate_symbols([n])
                expected = create_function_symbols(x, function=f)

                acctual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, acctual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2])
                expected = create_function_symbols(x, function=f)

                acctual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, acctual)

        # subtest 4
        # it returns [g_{0}, g_{1}, g_{2}, ...] when [a*h, b*h, c*h, ...] and 'g' are passed.
        num = random_int(2, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                acctual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, acctual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                acctual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, acctual)

        # subtest 5
        # it returns [f_{n}] when [n*h] and same_subscripts_as_stencil=True are passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                x = create_coordinate_symbols([n])
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                acctual = (sp.symbols(f + "_{%d}" % n),)
                self.assertEqual(expected, acctual)

                # staggered grid case
                stencil = [n / 2] if n != 0 else [n]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                str = f + "_{%2.1f}" % (n / 2) if n != 0 else f + "_{%d}" % (n)
                acctual = (sp.symbols(str),)
                self.assertEqual(expected, acctual)

        # subtest 6
        # it returns [f_{a}, f_{b}, f_{c}, ...] when [a*h, b*h, c*h, ...] and same_subscripts_as_stencil=True are passed.
        num = random_int(2, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                subscript = [
                    "_{%d}" % i if isinstance(i, int) else "_{%2.1f}" % i
                    for i in stencil
                ]
                str = "".join([f + s + " " for s in subscript])
                acctual = sp.symbols(str)
                self.assertEqual(expected, acctual)

                # staggered grid case
                stencil = [i / 2 if i != 0 else i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                subscript = [
                    "_{%d}" % i if isinstance(i, int) else "_{%2.1f}" % i
                    for i in stencil
                ]
                str = "".join([f + s + " " for s in subscript])
                acctual = sp.symbols(str)
                self.assertEqual(expected, acctual)

        # subtest 7
        # it returns [g_{n}] when [n*h], 'g', and same_subscripts_as_stencil=True are passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                x = create_coordinate_symbols([n])
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                acctual = (sp.symbols(f + "_{%d}" % n),)
                self.assertEqual(expected, acctual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2 if n != 0 else n])
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                str = f + "_{%2.1f}" % (n / 2) if n != 0 else f + "_{%d}" % (n)
                acctual = (sp.symbols(str),)
                self.assertEqual(expected, acctual)

        # subtest 8
        # it returns [g_{a}, g_{b}, g_{c}, ...] when [a*h, b*h, c*h, ...], 'g', and same_subscripts_as_stencil=True are passed.
        num = random_int(2, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                subscript = [
                    "_{%d}" % i if isinstance(i, int) else "_{%2.1f}" % i
                    for i in stencil
                ]
                str = "".join([f + s + " " for s in subscript])
                acctual = sp.symbols(str)
                self.assertEqual(expected, acctual)

                # staggered grid case
                stencil = [i / 2 if i != 0 else i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                subscript = [
                    "_{%d}" % i if isinstance(i, int) else "_{%2.1f}" % i
                    for i in stencil
                ]
                str = "".join([f + s + " " for s in subscript])
                acctual = sp.symbols(str)
                self.assertEqual(expected, acctual)

    def test_simplify_coefficients(self):
        """test suite for utils.simplify_coefficients.
        1. it returns integer list when sympy integer list is passed.
        2. it returns integer list when list of sympy symbols is passed.
        3. it returns list of sympy rational when sympy number list is passed.
        4. it returns list of sympy rational when list of sympy symbols is passed.
        5. it returns tuple of integer list and a integer when sympy integer list and as_numr_denom=True are passed.
        6. it returns tuple of integer list and a integer when lit of sympy symbols and as_numr_denom=True are passed.
        7.
        """

        # subtest 1
        # it returns sympy integers list when sympy integers list is passed.
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                numr = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                expected = [n / 1 for n in numr]

                acctual = simplify_coefficients(numr)

                self.assertEqual(expected, acctual)

        # subtest 2
        # it returns integer list when list of sympy symbols is passed.
        h = sp.symbols(random_string(5))
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                numr = [
                    random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    for _ in range(len_)
                ]
                coef = [n * h for n in numr]
                expected = [n / 1 for n in numr]

                acctual = simplify_coefficients(coef)

                self.assertEqual(expected, acctual)

        # subtest 3
        # it returns list of sympy rational when sympy number list is passed.
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                denom = [sp.Number(random.randint(1, 1000)) for _ in range(len_)]
                numr = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                coef = [numr[i] / denom[i] for i in range(len_)]
                expected = [sp.Rational(numr[i], denom[i]) for i in range(len_)]

                acctual = simplify_coefficients(coef)

                self.assertEqual(expected, acctual)

        # subtest 4
        # it returns list of sympy rational when list of sympy symbols is passed.
        h = sp.symbols(random_string(5))
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                denom = [sp.Number(random.randint(1, 1000)) for _ in range(len_)]
                numr = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                coef = [numr[i] / denom[i] * h for i in range(len_)]
                expected = [sp.Rational(numr[i], denom[i]) for i in range(len_)]

                acctual = simplify_coefficients(coef)

                self.assertEqual(expected, acctual)

        # subtest 5
        # it returns sympy integers list when sympy integers list is passed.
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                numr = [
                    sp.Number(
                        random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    )
                    for _ in range(len_)
                ]
                expected = ([n / 1 for n in numr], 1)

                acctual = simplify_coefficients(numr, as_numr_denom=True)

                self.assertEqual(expected, acctual)

        # subtest 6
        # it returns tuple of integer list and a integer when lit of sympy symbols and as_numr_denom=True are passed.
        h = sp.symbols(random_string(5))
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                numr = [
                    random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    for _ in range(len_)
                ]

                coef = [n * h for n in numr]
                expected = ([n / 1 for n in numr], 1)

                acctual = simplify_coefficients(coef, as_numr_denom=True)

                self.assertEqual(expected, acctual)

        # subtest 7
        for coef in [[sp.Number(3 / 4), sp.Number(2 / 3)]]:
            with self.subTest(coef):
                expected = (
                    [9, 8],
                    12,
                )

                acctual = simplify_coefficients(coef, as_numr_denom=True)

                self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
