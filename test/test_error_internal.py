"""Tests for distos.error.internal
"""
import sys

sys.path.insert(1, "..")

import unittest

from dictos.error.internal import UnexpectedDenominatorError


class ErrorInternalTest(unittest.TestCase):
    @unittest.expectedFailure
    def test_error_internal_UnexpectedDenominatorError_list(self):
        """
        test suite for error.internal.UnexpectedDenominatorError.
        """

        denom = [1]
        with self.subTest(denom):
            if (type(denom) is not list and type(denom) is not tuple) or (
                (type(denom) is list or type(denom) is tuple) and len(denom) > 1
            ):
                raise UnexpectedDenominatorError(denom)
            else:
                self.assertFalse(
                    (type(denom) is list or type(denom) is tuple) and len(denom) > 1
                )

        denom = [1, 2, 3]
        with self.subTest(denom):
            if (type(denom) is not list and type(denom) is not tuple) or (
                (type(denom) is list or type(denom) is tuple) and len(denom) > 1
            ):
                raise UnexpectedDenominatorError(denom)
            else:
                self.assertFalse(
                    (type(denom) is list or type(denom) is tuple) and len(denom) > 1
                )

    @unittest.expectedFailure
    def test_error_internal_UnexpectedDenominatorError_tuple(self):
        """
        test suite for error.internal.UnexpectedDenominatorError.
        """

        denom = (1,)
        with self.subTest(denom):
            if (type(denom) is not list and type(denom) is not tuple) or (
                (type(denom) is list or type(denom) is tuple) and len(denom) > 1
            ):
                raise UnexpectedDenominatorError(denom)
            else:
                self.assertFalse(
                    (type(denom) is list or type(denom) is tuple) and len(denom) > 1
                )

        denom = (1, 2, 3)
        with self.subTest(denom):
            if (type(denom) is not list and type(denom) is not tuple) or (
                (type(denom) is list or type(denom) is tuple) and len(denom) > 1
            ):
                raise UnexpectedDenominatorError(denom)
            else:
                self.assertFalse(
                    (type(denom) is list or type(denom) is tuple) and len(denom) > 1
                )

    @unittest.expectedFailure
    def test_error_internal_UnexpectedDenominatorError_int(self):
        """
        test suite for error.internal.UnexpectedDenominatorError.
        """

        denom = 1
        with self.subTest(denom):
            if (type(denom) is not list and type(denom) is not tuple) or (
                (type(denom) is list or type(denom) is tuple) and len(denom) > 1
            ):
                raise UnexpectedDenominatorError(denom)
            else:
                self.assertFalse(
                    (type(denom) is list or type(denom) is tuple) and len(denom) > 1
                )


if __name__ == "__main__":
    unittest.main()
