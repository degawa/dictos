"""Tests for test.gen
"""

import sys

sys.path.insert(1, "..")

import unittest
import random

from test.utilities.gen import random_int, random_string

_MAX_LENGTH = 100


class GenTest(unittest.TestCase):
    def test_random_int(self):
        """test suite for test.random_int."""

        min = random.randint(-_MAX_LENGTH, 0)
        max = random.randint(0, _MAX_LENGTH)

        num = random_int(min, max)
        ordered_list = list(range(min, max + 1))

        with self.subTest("random_int returns unordered list"):
            self.assertTrue(len(num) == len(ordered_list) and num != ordered_list)

        with self.subTest(
            "random_int returns a list in which the same number does not appear twice."
        ):
            self.assertTrue(len(num) == len(set(num)))

        exclude = [0]
        num = random_int(min, max, exclude)
        with self.subTest(
            "random_int returns a list not including numbers specified by `exclude`"
        ):
            for i in exclude:
                self.assertTrue(i not in num)

        with self.subTest(
            "random_int returns unordered list not including numbers specified by `exclude`"
        ):
            for i in exclude:
                ordered_list.remove(i)

            self.assertTrue(len(num) == len(ordered_list) and num != ordered_list)

        with self.subTest(
            "random_int returns an excluded list in which the same number does not appear twice."
        ):
            self.assertTrue(len(num) == len(set(num)))

    def test_random_string(self):
        """test suite for test.random_string."""

        strlen = random.randint(0, _MAX_LENGTH)

        string = random_string(strlen)

        with self.subTest("random_string returns str"):
            self.assertTrue(type(string) is str)

        with self.subTest(f"random_string returns str that has length {strlen}"):
            self.assertTrue(len(string) == strlen)


if __name__ == "__main__":
    unittest.main()
