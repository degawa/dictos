"""Tests for distos.utils
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.utils import (
    DEFAULT_INTERVAL_SYMBOL_STR,
    create_set_of_coordinate_symbols_from_stencil,
    create_set_of_function_symbols_at_coordinate,
    simplify_coefficients,
    dotproduct,
    div,
)

_stencil_half_width = 20  # up to 20th order accuracy


class UtilsTest(unittest.TestCase):
    def test_create_set_of_coordinate_symbols_from_stencil(self):
        """
        test suite for utils.create_set_of_coordinate_symbols_from_stencil.
        it returns [nh] when [n] is passed.
        """

        num = list(range(-_stencil_half_width, 0))
        num += list(range(1, _stencil_half_width + 1))
        random.shuffle(num)
        for n in num:
            with self.subTest(n):
                stencil = [n]
                expected = create_set_of_coordinate_symbols_from_stencil(stencil)

                h = sp.symbols(DEFAULT_INTERVAL_SYMBOL_STR)
                acctual = [n * h]
                self.assertEqual(expected, acctual)

        num = list(range(1, _stencil_half_width + 1))
        random.shuffle(num)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                expected = create_set_of_coordinate_symbols_from_stencil(stencil)

                h = sp.symbols(DEFAULT_INTERVAL_SYMBOL_STR)
                acctual = [i * h for i in range(n)]
                self.assertEqual(expected, acctual)

        # uncomment test after implement error handling
        # with self.subTest():
        #     with self.assertRaises(Exception):
        #         create_set_of_coordinate_symbols_from_stencil([])


if __name__ == "__main__":
    unittest.main()
