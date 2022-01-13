"""Tests for distos.utils
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random
import string

from dictos.utils import (
    DEFAULT_INTERVAL,
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
    """generate random list of integers from from_ to to_.
    Integers listed in exclude_ are excluded

    Args:
        min (integer): minimum value of integer in the list.
        max (integer): maximum value of integer in the list.
        exclude (list of integer): excluded integers from the list.
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


if __name__ == "__main__":
    unittest.main()
