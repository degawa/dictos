"""Tests for distos.interplation
"""
import sys

sys.path.insert(1, "..")

import unittest
import sympy as sp

from dictos.interpolation import equation, coefficients, truncation_error


class UtilsTest(unittest.TestCase):
    def test_equation(self):
        """
        test suite for interplation.equation.
        """

        f_0 = sp.symbols("f_0")
        f_1 = sp.symbols("f_1")
        f_2 = sp.symbols("f_2")
        f_3 = sp.symbols("f_3")
        f_4 = sp.symbols("f_4")
        f_5 = sp.symbols("f_5")
        f_6 = sp.symbols("f_6")
        f_7 = sp.symbols("f_7")
        f_8 = sp.symbols("f_8")
        f_9 = sp.symbols("f_9")
        f_10 = sp.symbols("f_10")

        expected = [
            0,
            f_0 / 2 + f_1 / 2,
            -f_0 / 6 + 2 * f_1 / 3 + 2 * f_2 / 3 - f_3 / 6,
            f_0 / 20
            - 3 * f_1 / 10
            + 3 * f_2 / 4
            + 3 * f_3 / 4
            - 3 * f_4 / 10
            + f_5 / 20,
            -f_0 / 70
            + 4 * f_1 / 35
            - 2 * f_2 / 5
            + 4 * f_3 / 5
            + 4 * f_4 / 5
            - 2 * f_5 / 5
            + 4 * f_6 / 35
            - f_7 / 70,
            f_0 / 252
            - 5 * f_1 / 126
            + 5 * f_2 / 28
            - 10 * f_3 / 21
            + 5 * f_4 / 6
            + 5 * f_5 / 6
            - 10 * f_6 / 21
            + 5 * f_7 / 28
            - 5 * f_8 / 126
            + f_9 / 252,
        ]
        for half_width in range(1, 6):
            with self.subTest("%d-point central interpolation" % (half_width * 2)):
                stencil = [i for i in range(-half_width, half_width + 1)]
                stencil.remove(0)
                actual = equation(stencil)
                self.assertEqual(expected[half_width], actual)

        expected = [
            0,
            0,
            2 * f_0 - f_1,
            3 * f_0 - 3 * f_1 + f_2,
            4 * f_0 - 6 * f_1 + 4 * f_2 - f_3,
            5 * f_0 - 10 * f_1 + 10 * f_2 - 5 * f_3 + f_4,
            6 * f_0 - 15 * f_1 + 20 * f_2 - 15 * f_3 + 6 * f_4 - f_5,
            7 * f_0 - 21 * f_1 + 35 * f_2 - 35 * f_3 + 21 * f_4 - 7 * f_5 + f_6,
            8 * f_0
            - 28 * f_1
            + 56 * f_2
            - 70 * f_3
            + 56 * f_4
            - 28 * f_5
            + 8 * f_6
            - f_7,
            9 * f_0
            - 36 * f_1
            + 84 * f_2
            - 126 * f_3
            + 126 * f_4
            - 84 * f_5
            + 36 * f_6
            - 9 * f_7
            + f_8,
        ]
        for width in range(2, 10):
            with self.subTest("%d-point forward interpolation" % (width * 2)):
                stencil = [i for i in range(1, width + 1)]
                actual = equation(stencil)
                self.assertEqual(expected[width], actual)



if __name__ == "__main__":
    unittest.main()
