"""Tests for distos.error.stencil
"""
import sys

sys.path.insert(1, "..")

import unittest

from dictos.utils import has_zero
from dictos.error.stencil import ContainsZeroError, DuplicatedPointError


class ErrorStencilTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_stencil_ContainsZeroError(self):
        """
        test suite for error.stencil.ContainsZeroError.
        """

        stencil = [0, 1, 2]
        if has_zero(stencil):
            raise ContainsZeroError


if __name__ == "__main__":
    unittest.main()
