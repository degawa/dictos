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
    create_set_of_function_symbols_at_coordinate,
    simplify_coefficients,
    dotproduct,
    div,
)

_stencil_half_width = 20  # up to 20th order accuracy


def random_string(len):
    """generate n-length random string

    Args:
        len (integer): length of string
    """

    return "".join(random.choices(string.ascii_letters + string.digits, k=len))


class UtilsTest(unittest.TestCase):
    def test_create_set_of_coordinate_symbols_from_stencil(self):
        """
        test suite for utils.create_set_of_coordinate_symbols_from_stencil.
        - it returns [nh] when [n] is passed.
        - it returns [ah bh ch ...] when [a b c ...] is passed.
        - it returns [adx bdx cdx ...] when [a b c ...] and 'dx' are passed.
        """

        num = list(range(-_stencil_half_width, 0)) + list(
            range(1, _stencil_half_width + 1)
        )
        random.shuffle(num)
        for n in num:
            with self.subTest(n):
                stencil = [n]
                expected = create_coordinate_symbols(stencil)

                h = sp.symbols(DEFAULT_INTERVAL)
                acctual = [n * h]
                self.assertEqual(expected, acctual)

        num = list(range(1, _stencil_half_width + 1))
        random.shuffle(num)
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

        num = list(range(1, _stencil_half_width + 1))
        random.shuffle(num)
        for n in num:
            with self.subTest(n):
                interval = random_string(random.randint(1, 10))
                h = sp.symbols(interval)

                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil, interval=interval)

                acctual = [i * h for i in range(n)]
                self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
