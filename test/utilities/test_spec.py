"""Tests for distos.utilities.spec
"""

import sys

sys.path.insert(1, "..")

import unittest
import random

from dictos.utilities.spec import (
    has_zero,
    has_duplicated_points,
    narrower_than_minimum_width,
    is_not_natural_number,
    is_not_positive_integer,
    is_not_assumed_length,
    are_different_length,
    _MINIMUM_STENCIL_WIDTH,
)
from test.utilities.gen import random_int


class SpecTest(unittest.TestCase):
    def test_has_zero(self):
        """test suite for spec.has_zero.
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
        """test suite for spec.has_dupulicated_points.
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

    def test_narrower_than_minimum_width(self):
        """test suite for spec.narrower_than_minimum_width."""

        for c in random_int(-100, 100):
            with self.subTest(c):
                stencil = list(range(c))
                expected = True if c < _MINIMUM_STENCIL_WIDTH else False
                actual = narrower_than_minimum_width(stencil)
                self.assertEqual(expected, actual)

    def test_is_not_natrual_number(self):
        """test suite for spec.is_not_natrual_number."""

        for c in random_int(-100, 100):
            with self.subTest(c):
                expected = True if c <= 0 else False
                actual = is_not_natural_number(c)
                self.assertEqual(expected, actual)

            with self.subTest(c):
                expected = True if c < 1 else False
                actual = is_not_natural_number(c)
                self.assertEqual(expected, actual)

    def test_is_not_positive_integer(self):
        """test suite for spec.is_not_positive_integer."""

        for c in random_int(-100, 100):
            with self.subTest(c):
                expected = True if c <= 0 else False
                actual = is_not_natural_number(c)
                self.assertEqual(expected, actual)

            with self.subTest(c):
                expected = True if c < 0 else False
                actual = is_not_positive_integer(c, include_zero=True)
                self.assertEqual(expected, actual)

    def test_is_not_assumed_length(self):
        """test suite for spec.is_not_assumed_length."""

        for c in random_int(1, 100):
            with self.subTest(c):
                stencil = list(range(c))
                self.assertTrue(is_not_assumed_length(stencil, c + 1))

            with self.subTest(c):
                stencil = list(range(c))
                self.assertFalse(is_not_assumed_length(stencil, c))

    def test_are_different_length(self):
        """test suite for spec.are_different_length."""

        for c in random_int(1, 100):
            with self.subTest(c):
                list1 = list(range(c))
                list2 = list(range(c))
                self.assertFalse(are_different_length(list1, list2))

            with self.subTest(c):
                list1 = list(range(c + 1))
                list2 = list(range(c))
                self.assertTrue(are_different_length(list1, list2))

                list1 = list(range(c))
                list2 = list(range(c + 1))
                self.assertTrue(are_different_length(list1, list2))


if __name__ == "__main__":
    unittest.main()
