"""Tests for distos.error.linear_algebra
"""
import sys

sys.path.insert(1, "..")

import unittest

from dictos.error.linear_algebra import InconsistentDataSetError


class ErrorLinearAlgebraTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_linear_algebra_InconsistentDataSetError(self):
        """
        test suite for error.linear_algebra.InconsistentDataSetError.
        """

        vec1 = [1, 2, 3]
        vec2 = [1, 2, 3, 4]
        with self.subTest(vec1, vec2):
            if len(vec1) != len(vec2):
                raise InconsistentDataSetError(vec1, vec2)


if __name__ == "__main__":
    unittest.main()
