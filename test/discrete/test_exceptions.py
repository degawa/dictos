"""Tests for distos.discrete.exceptions
"""

import sys

sys.path.insert(1, "..")

import unittest

from dictos.utilities.spec import has_zero, has_duplicated_points
from dictos.discrete.exceptions import (
    ContainsZeroError,
    DuplicatedPointError,
    TooNarrowError,
)


class ErrorStencilTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_stencil_ContainsZeroError(self):
        """
        test suite for discrete.exceptions.ContainsZeroError.
        """

        stencil = [0, 1, 2]
        if has_zero(stencil):
            raise ContainsZeroError

    @unittest.expectedFailure
    def test_error_stencil_DuplicatedPointError(self):
        """
        test suite for discrete.exceptions.DuplicatedPointError.
        """

        with self.subTest():
            stencil = [0, 1, 2]
            if has_duplicated_points(stencil):
                raise DuplicatedPointError(stencil)

        with self.subTest():
            stencil = [1, 2, 1]
            if has_duplicated_points(stencil):
                raise DuplicatedPointError(stencil)

    @unittest.expectedFailure
    def test_error_stencil_TooNarrowError(self):
        """
        test suite for discrete.exceptions.TooNarrowError.
        """

        with self.subTest():
            stencil = [0, 1]
            if len(stencil) <= 1:
                raise TooNarrowError(stencil)

        with self.subTest():
            stencil = [1]
            if len(stencil) <= 1:
                raise TooNarrowError(stencil)


if __name__ == "__main__":
    unittest.main()
