"""Tests for distos.spec
"""
import sys

sys.path.insert(1, "..")

import unittest
import random

from dictos.spec import (
    has_zero,
    has_duplicated_points,
)
from gen import random_int


class SpecTest(unittest.TestCase):
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
