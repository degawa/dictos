"""Tests for distos.stencil
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp
import random

from dictos.spec import DEFAULT_INTERVAL, DEFAULT_FUNCTION
from dictos.stencil import create_coordinate_symbols, create_function_symbols
from gen import random_string, random_int, STENCIL_HALF_WIDTH, MAX_SYMBOL_LENGTH


class StencilTest(unittest.TestCase):
    def test_create_coordinate_symbols(self):
        """
        test suite for stencil.create_coordinate_symbols.
        1. it returns [a*h b*h c*h ...] when [a b c ...] is passed.
        2. it returns [a*dx b*dx c*dx ...] when [a b c ...] and 'dx' are passed.
        3. it raise error when empty list is passed.
        4. it raise error when at least a number in the stencil appears more than once.
        """

        # subtest 1
        # it returns [a*h b*h c*h ...] when [a b c ...] is passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil)

                h = sp.symbols(DEFAULT_INTERVAL)
                actual = [i * h for i in range(n)]
                self.assertEqual(expected, actual)

        # subtest 2
        # it returns [a*dx b*dx c*dx ...] when [a b c ...] and 'dx' are passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                interval = random_string(random.randint(1, MAX_SYMBOL_LENGTH))
                h = sp.symbols(interval)

                stencil = [i for i in range(n)]
                expected = create_coordinate_symbols(stencil, interval=interval)

                actual = [i * h for i in range(n)]
                self.assertEqual(expected, actual)

        # subtest 3 & 4 are tested in `test_error_stencil` moudle

    def test_create_function_symbols(self):
        """
        test suite for stencil.create_function_symbols.
        1. it returns [f_{0}, f_{1}, f_{2}, ...] when [a*h, b*h, c*h, ...] is passed.
        2. it returns [g_{0}, g_{1}, g_{2}, ...] when [a*h, b*h, c*h, ...] and 'g' are passed.
        3. it returns [f_{a}, f_{b}, f_{c}, ...] when [a*h, b*h, c*h, ...] and same_subscripts_as_stencil=True are passed.
        4. it returns [g_{a}, g_{b}, g_{c}, ...] when [a*h, b*h, c*h, ...], 'g', and same_subscripts_as_stencil=True are passed.
        """

        # subtest 1
        # it returns [f_{0}, f_{1}, f_{2}, ...] when [a*h, b*h, c*h, ...] is passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                actual = sp.symbols(f + f"_0:{n}")
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x)

                f = DEFAULT_FUNCTION
                actual = sp.symbols(f + f"_0:{n}")
                self.assertEqual(expected, actual)

        # subtest 2
        # it returns [g_{0}, g_{1}, g_{2}, ...] when [a*h, b*h, c*h, ...] and 'g' are passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, MAX_SYMBOL_LENGTH))

                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                actual = sp.symbols(f + f"_0:{n}")
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, function=f)

                actual = sp.symbols(f + f"_0:{n}")
                self.assertEqual(expected, actual)

        # subtest 3
        # it returns [f_{a}, f_{b}, f_{c}, ...] when [a*h, b*h, c*h, ...] and same_subscripts_as_stencil=True are passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                f = DEFAULT_FUNCTION
                subscript = [
                    f"{i}" if isinstance(i, int) else f"{i:2.1f}" for i in stencil
                ]
                str = "".join([f + "_{" + s + "}" + " " for s in subscript])
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 if i != 0 else i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(x, same_subscripts_as_stencil=True)

                subscript = [
                    f"{i}" if isinstance(i, int) else f"{i:2.1f}" for i in stencil
                ]
                str = "".join([f + "_{" + s + "}" + " " for s in subscript])
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

        # subtest 4
        # it returns [g_{a}, g_{b}, g_{c}, ...] when [a*h, b*h, c*h, ...], 'g', and same_subscripts_as_stencil=True are passed.
        num = random_int(2, STENCIL_HALF_WIDTH)
        for n in num:
            with self.subTest(n):
                f = random_string(random.randint(1, MAX_SYMBOL_LENGTH))

                stencil = [i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                subscript = [
                    f"{i}" if isinstance(i, int) else f"{i:2.1f}" for i in stencil
                ]
                str = "".join([f + "_{" + s + "}" + " " for s in subscript])
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)

                # staggered grid case
                stencil = [i / 2 if i != 0 else i for i in range(n)]
                x = create_coordinate_symbols(stencil)
                expected = create_function_symbols(
                    x, function=f, same_subscripts_as_stencil=True
                )

                subscript = [
                    f"{i}" if isinstance(i, int) else f"{i:2.1f}" for i in stencil
                ]
                str = "".join([f + "_{" + s + "}" + " " for s in subscript])
                actual = sp.symbols(str)
                self.assertEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()