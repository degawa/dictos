"""Tests for distos.linalg.exceptions
"""

import sys

sys.path.insert(1, "..")

import unittest

from dictos.linalg.exceptions import InconsistentDataSetError


class ErrorLinearAlgebraTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_linear_algebra_InconsistentDataSetError(self):
        """
        test suite for linalg.exceptions.InconsistentDataSetError.
        """

        vec1 = [1, 2, 3]
        vec2 = [1, 2, 3, 4]
        with self.subTest(expected=vec1, actual=vec2):
            if len(vec1) != len(vec2):
                raise InconsistentDataSetError(vec1, vec2)


if __name__ == "__main__":
    unittest.main()
