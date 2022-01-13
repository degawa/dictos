import sympy as sp
import numpy as np
import fractions as fr


DEFAULT_INDEPENDENT_VARIABLE = "x"  # str for independent variable symbol
DEFAULT_INTERVAL = "h"  # str for interval symbol
DEFAULT_FUNCTION = "f"  # str for function symbol


def create_coordinate_symbols(stencil, interval=DEFAULT_INTERVAL):
    """
    create set of coordinate symbols from stencil.
    input a list of numbers like `[-1, 0, 1]` as stencil,
    this returns a list of coordinates like `[-h, 0, h]`.

    Args:
        stencil (list of int or float): stencil on regular or
            staggered grid. It is not allowed that a number
            in the list appears more than once.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL_SYMBOL_STR.

    Returns:
        list of sympy symbols: list of coordinates
            corresponding to the stencil.
    """

    # TODO: raise error when len(stencil)==0
    # TODO: raise error when at least a number in the stencil appears more than once.

    return [stencil[i] * sp.symbols(interval) for i in range(len(stencil))]


def create_function_symbols(
    x_set,
    function=DEFAULT_FUNCTION,
    same_subscripts_as_stencil=False,
):
    """
    create set of function symbols at coordinates from x_set.
    input a list of sympy symbols like `[-h, 0, h]` as x_set,
    this returns a list of functions at x_set like `[f_0, f_1, f_2]`.

    Args:
        x_set (list of sympy symbols): coordinates on regular or
            staggered grid corresponding to the stencil.
        function (str, optional): a function symbol like `f`.
            Defaults to DEFAULT_FUNCTION.
        same_subscripts_as_stencil (bool, optional): flag
            to make function subscripts the same as the stencil.
            Defaults to False, the subscripts start from 0.

    Returns:
        list of sympy symbols: list of functions at passed coordinates
    """
    if same_subscripts_as_stencil:
        # make function subscripts the same as the stencil

        stencil = [x if x.is_number else sp.poly(x).coeffs()[0] for x in x_set]
        # extruct numbers from list of coordinates.
        # coordinate consists of a number and a symbol,
        # such as -2*h and 1.5*h, extract the number as coefficint.
        # coordinate is 0, that is a number, use 0 as the stencil

        subscript = ["_{%d}" % i if i.is_integer else "_{%2.1f}" % i for i in stencil]
        str = "".join([function + s + " " for s in subscript])
        # make string "f_{-1} f_{-0.5} ..." to pass `sympy.symbols`.
        # tail space is ignored in `sympy.symbols`

        f_set = sp.symbols(str)
    else:
        f_set = sp.symbols(function + "_0:{:d}".format(len(x_set)))
        # make a list of sympy symbols from string "f_{0} f_{1} ...".

    return f_set


def simplify_coefficients(coef, as_numr_denom):
    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]

    coef_rational = [
        sp.Rational(fr.Fraction(str(c)).limit_denominator(100000)) for c in coef_num
    ]

    denom = [c.q for c in coef_rational]
    denom_lcm = np.lcm.reduce(np.array(denom))
    numr = [c * denom_lcm for c in coef_rational]

    if as_numr_denom:
        return numr, denom_lcm
    else:
        return [n / denom_lcm for n in numr]


def dotproduct(numr, f_set):
    idx_start = len(f_set) - 1
    idx_end = -1
    eq = sp.Mul(numr[idx_start], f_set[idx_start], evaluate=False)
    for i in range(idx_start - 1, idx_end, -1):
        eq = sp.Add(eq, sp.Mul(numr[i], f_set[i], evaluate=False), evaluate=False)

    return eq


def div(eq, denom):
    return sp.Mul(eq, 1 / denom)
