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


class UtilsTest(unittest.TestCase):
    def test_create_set_of_coordinate_symbols_from_stencil(self):
        """
        test suite for utils.create_set_of_coordinate_symbols_from_stencil.
        it returns [nh] when [n] is passed.
        """

        num = range(-20, 21, 1)
        random.shuffle(num)
        for n in num:
            with self.subTest(n):
                stencil = [n]
                expected = create_set_of_coordinate_symbols_from_stencil(stencil)

                h = sp.symbols(DEFAULT_INTERVAL_SYMBOL_STR)
                acctual = [n * h]
                self.assertEqual(expected, acctual)


if __name__ == "__main__":
    unittest.main()
