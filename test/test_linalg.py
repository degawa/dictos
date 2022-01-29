"""Tests for distos.linalg
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.linalg import dot_product, div
from dictos.error.linear_algebra import InconsistentDataSetError
from gen import random_int, random_string


class LinalgTest(unittest.TestCase):
    def test_linalg_dot_product(self):
        """test suite for linalg.dot_product"""

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            vec1 = random_int(-n, n)
            vec2 = random_int(-n, n)
            with self.subTest(
                f"dot_product of 2 vectors with {len(vec1)} elements with evaluate True"
            ):
                expected = sum([vec1[i] * vec2[i] for i in range(len(vec1))])
                actual = dot_product(vec1, vec2, evaluate=True)
                self.assertEqual(expected, actual)

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            vec1 = [sp.symbols(random_string(2)) for _ in range(num_max)]
            vec2 = [sp.symbols(random_string(2)) for _ in range(num_max)]
            with self.subTest(
                f"dot_product of 2 vectors with {len(vec1)} elements with evaluate True"
            ):
                expected = sum([vec1[i] * vec2[i] for i in range(len(vec1))])
                actual = dot_product(vec1, vec2, evaluate=True)
                self.assertEqual(expected, actual)

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            vec1 = [sp.symbols(random_string(2)) for _ in range(num_max)]
            vec2 = [sp.symbols(random_string(2)) for _ in range(num_max)]
            with self.subTest(
                f"dot_product of 2 vectors with {len(vec1)} elements with evaluate False"
            ):
                expected = sum([vec1[i] * vec2[i] for i in range(len(vec1))])
                actual = dot_product(vec1, vec2, evaluate=False)
                self.assertNotEqual(expected, actual)

        num_max = 3
        # for n in range(1, num_max):
        vec1 = [sp.symbols(random_string(2)) for _ in range(num_max)]
        vec2 = [sp.symbols(random_string(2)) for _ in range(num_max)]
        with self.subTest(
            f"dot_product of 2 vectors with {len(vec1)} elements with evaluate False"
        ):
            begin_ = len(vec2) - 1
            end_ = -1
            step_ = -1
            expected = sp.Mul(vec1[begin_], vec2[begin_], evaluate=False)
            for i in range(begin_ + step_, end_, step_):
                expected = sp.Add(
                    expected, sp.Mul(vec1[i], vec2[i], expected=False), evaluate=False
                )
            actual = dot_product(vec1, vec2, evaluate=False)

            ac_str = str(actual)
            ex_str = str(expected)
            self.assertTrue(
                sp.simplify(expected) == sp.simplify(actual) and ac_str == ex_str
            )
            # I couldn't find a way to get True when evaluate=False without simplify.
            # string converted from sympy expr looks same as original expr.

    def test_linalg_div(self):
        """test suite for linalg.dot_product"""

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            numer = random_int(-n, n)
            denom = random.randint(1, 100)
            with self.subTest(f"dev of {len(numer)} elements numerator and {denom}"):
                expected = sum([n for n in numer]) / denom
                actual = div(sum([n for n in numer]), denom)
                self.assertEqual(expected, actual)

        num_max = random.randint(1, 10)
        for n in range(1, num_max):
            numer = [sp.symbols(random_string(2)) for _ in range(num_max)]
            denom = sp.symbols(random_string(2))
            with self.subTest(f"dev of {len(numer)} elements numerator and {denom}"):
                expected = sum([n for n in numer]) / denom
                actual = div(sum([n for n in numer]), denom)
                self.assertEqual(expected, actual)

    def test_linalg_exception(self):
        """test suite for exception in linalg"""
        numer = [sp.core.Number(i) for i in range(4)]
        f_set = sp.symbols(f"f_0:{(len(numer) + 1)}")
        with self.subTest("dot_product with inconsistent numer and f_set"):
            with self.assertRaises(InconsistentDataSetError):
                dot_product(numer, f_set)


if __name__ == "__main__":
    unittest.main()
