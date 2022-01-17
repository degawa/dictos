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
    has_zero,
    has_duplicated_points,
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
                actual = [n * h]
                self.assertEqual(expected, actual)

        # subtest 2
        # it returns [a*h b*h c*h ...] when [a b c ...] is passed.
        num = random_int(1, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil)

                h = sp.symbols(DEFAULT_INTERVAL)
                actual = [i * h for i in range(n)]
                self.assertEqual(expected, actual)

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

                actual = [i * h for i in range(n)]
                self.assertEqual(expected, actual)

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
                actual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, actual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2])
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                actual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, actual)

        # subtest 2
        # it returns [f_{0}, f_{1}, f_{2}, ...] when [a*h, b*h, c*h, ...] is passed.
        num = random_int(1, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                actual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                actual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, actual)

        # subtest 3
        # it returns [g_{0}] when [n*h] and 'g' are passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                x = create_coordinate_symbols([n])
                expected = create_function_symbols(x, function=f)

                actual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, actual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2])
                expected = create_function_symbols(x, function=f)

                actual = (sp.symbols(f + "_0"),)
                self.assertEqual(expected, actual)

        # subtest 4
        # it returns [g_{0}, g_{1}, g_{2}, ...] when [a*h, b*h, c*h, ...] and 'g' are passed.
        num = random_int(2, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, _max_symbol_length))

                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                actual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                actual = sp.symbols((f + "_0:{:d}").format(n))
                self.assertEqual(expected, actual)

        # subtest 5
        # it returns [f_{n}] when [n*h] and same_subscripts_as_stencil=True are passed.
        num = random_int(-_stencil_half_width, _stencil_half_width)
        for n in num:
            with self.subTest(n):
                x = create_coordinate_symbols([n])
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                actual = (sp.symbols(f + "_{%d}" % n),)
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [n / 2] if n != 0 else [n]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                str = f + "_{%2.1f}" % (n / 2) if n != 0 else f + "_{%d}" % (n)
                actual = (sp.symbols(str),)
                self.assertEqual(expected, actual)

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
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 if i != 0 else i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                subscript = [
                    "_{%d}" % i if isinstance(i, int) else "_{%2.1f}" % i
                    for i in stencil
                ]
                str = "".join([f + s + " " for s in subscript])
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

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

                actual = (sp.symbols(f + "_{%d}" % n),)
                self.assertEqual(expected, actual)

                # staggered grid case
                x = create_coordinate_symbols([n / 2 if n != 0 else n])
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                str = f + "_{%2.1f}" % (n / 2) if n != 0 else f + "_{%d}" % (n)
                actual = (sp.symbols(str),)
                self.assertEqual(expected, actual)

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
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

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
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

    def test_simplify_coefficients(self):
        """test suite for utils.simplify_coefficients.
        1. it returns integer list when sympy integer list is passed.
        2. it returns integer list when list of sympy symbols is passed.
        3. it returns list of sympy rational when sympy number list is passed.
        4. it returns list of sympy rational when list of sympy symbols is passed.
        5. it returns tuple of integer list and a integer when sympy integer list and as_numr_denom=True are passed.
        6. it returns tuple of integer list and a integer when list of sympy symbols and as_numr_denom=True are passed.
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

                actual = simplify_coefficients(numr)

                self.assertEqual(expected, actual)

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

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

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

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

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

                actual = simplify_coefficients(coef)

                self.assertEqual(expected, actual)

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

                actual = simplify_coefficients(numr, as_numr_denom=True)

                self.assertEqual(expected, actual)

        # subtest 6
        # it returns tuple of integer list and a integer when list of sympy symbols and as_numr_denom=True are passed.
        h = sp.symbols(random_string(5))
        for len_ in random_int(1, _stencil_half_width * 2):
            with self.subTest(len_):
                numr = [
                    random.randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max)
                    for _ in range(len_)
                ]

                coef = [n * h for n in numr]
                expected = ([n / 1 for n in numr], 1)

                actual = simplify_coefficients(coef, as_numr_denom=True)

                self.assertEqual(expected, actual)

        # subtest 7
        for coef in [[sp.Number(3 / 4), sp.Number(2 / 3), sp.Number(4 / 5)]]:
            with self.subTest(coef):
                expected = (
                    [45, 40, 48],
                    60,
                )

                actual = simplify_coefficients(coef, as_numr_denom=True)

                self.assertEqual(expected, actual)

    def test_has_zero(self):
        """test suite for utils.has_zero.
        1. returns True if stencil has 0
        2. returns False if stencil doesn't have 0
        """

        # 1
        for c in range(100):
            stencil = random_int(-c, c)
            with self.subTest(stencil):
                expected = True
                actual = has_zero(stencil)
                self.assertEqual(expected, actual)

        # 2
        for c in range(100):
            stencil = random_int(-c, c, exclude=[0])
            with self.subTest(stencil):
                expected = False
                actual = has_zero(stencil)
                self.assertEqual(expected, actual)

    def test_has_duplicated_points(self):
        """test suite for utils.has_dupulicated_points.
        1. returns True if there is at least one dupulicated point.
        2. returns False if there is no dupulicated point.
        """

        # 1
        for c in range(100):
            stencil = random_int(-c, c)
            stencil.append(random.randint(-c, c))
            random.shuffle(stencil)
            with self.subTest(stencil):
                expected = True
                actual = has_duplicated_points(stencil)
                self.assertEqual(expected, actual)

        # 2
        for c in range(100):
            stencil = random_int(-c, c)
            with self.subTest(stencil):
                expected = False
                actual = has_duplicated_points(stencil)
                self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
