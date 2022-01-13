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
        tuple of sympy symbols: tuple of functions at passed coordinates.
            When same_subscripts_as_stencil is set,
            subscripts is enclused in curly braces like `f_{0}`.
    """
    if same_subscripts_as_stencil:
        # make function subscripts the same as the stencil

        stencil = [x if x.is_number else sp.poly(x).coeffs()[0] for x in x_set]
        # extract numbers from list of coordinates.
        # coordinate consists of a number and a symbol,
        # such as -2*h and 1.5*h, extract the number as coefficint.
        # coordinate is 0, that is a number, use 0 as the stencil

        subscript = ["_{%d}" % i if i.is_integer else "_{%2.1f}" % i for i in stencil]
        str = "".join([function + s + " " for s in subscript])
        # make string "f_{-1} f_{-0.5} ..." to pass `sympy.symbols`.
        # tail space is ignored in `sympy.symbols`

        f_set = sp.symbols(str)
        if type(f_set) == sp.core.symbol.Symbol:
            f_set = (f_set,)
        # make the return value's type a tuple.
        # sp.symbols(str) returns sp.core.symbol.Symbol
        # when str does not contain space
        # althought sp.symbols("_0:{:d}".format(n)) returns tuple
        # when n==1.
    else:
        f_set = sp.symbols(function + "_0:{:d}".format(len(x_set)))
        # make a tuple of sympy symbols from string.

    return f_set


def simplify_coefficients(coef, as_numr_denom=False):
    """
    simplify coefficients in floating-point number
    as ratioanl numbers.
    When as_numr_denom flag is set, this returns the numberator and
    denominator separately.

    Args:
        coef (list of sympy Mul and numbers): finite difference and
            interpolation coefficients to be simplified.
        as_numr_denom (bool, optional): flag to return the numerator
            and denominator separately.
            Defaults to False.

    Returns:
        list of sympy Rational: simplified coefficients.\\
        list of sympy numbers, numpy.int32:
            numerator of simplified coefficients,
            and denominator as the least common multiple
    """
    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]
    # extract numbers from list of coefficients like [1/h**2, ...].

    coef_rational = [
        sp.Rational(fr.Fraction(str(c)).limit_denominator(100000)) for c in coef_num
    ]
    # rationalize coefficients.
    # The argument of `fraction.Fraction` must be a string in this case.
    # `limit_denominator` is used to avoid represent a real number
    # including the error coused by the binary representation.
    # It is necessary to discuss the argument
    # when increasing the formal accuracy of discretization.
    # TODO: calculate the argument of `limit_denominator` from the smallest coefficient, that is the largest denominator

    denom = [c.q for c in coef_rational]
    denom_lcm = np.lcm.reduce(np.array(denom))
    # extract denomenator of each coefficient
    # and calculate the least common multiple.

    numr = [c * denom_lcm for c in coef_rational]
    # list of numerator divided by the least common multiple.

    if as_numr_denom:
        return numr, denom_lcm
    else:
        return [n / denom_lcm for n in numr]


def dot_product(numr, f_set, evaluate=False):
    """
    calculate dot product of two list
    containing sympy numbers and symbols without evaluation.

    Args:
        numr (list of sympy Mul and numbers):
            list of numerator of coefficients.
            The list length must be equal to it of f_set.
        f_set (list of sympy symbols): list of function symbols.
            The list length must be equal to it of numr.
        evaluate (bool, optional): flag to evaluate the result.
            Defaults to False.

    Returns:
        sympy Expr: dot product of the passed two lists.
    """

    # TODO: raise error when numr and f_set have different lengths

    begin_ = len(f_set) - 1  # exclude first term
    end_ = -1  # to generate numbers up to 0 using range()
    step_ = -1
    eq = sp.Mul(numr[begin_], f_set[begin_], evaluate=evaluate)
    for i in range(begin_ + step_, end_, step_):
        eq = sp.Add(eq, sp.Mul(numr[i], f_set[i], evaluate=evaluate), evaluate=evaluate)
    # there is no way to accumulate a list of sympy symbols
    # without evalulation.
    # So for-loop and sympy.Add and .Mul are used.
    # The for-loop is downstepped
    # so that the result are sorted when it is printed.
    # The result of the accumulation, `((-a*f_{-1} +b*f_{0}) + a*f_{1})`,
    # is printed as `a*f_{1} + b*f_{0} - a*f_{-1}`.
    # downstepping is used for sorting.

    return eq


def div(eq, denom):
    """
    calculate division of two sympy Expr.
    The result is always evaluated.

    Args:
        eq (sympy Expr): Algebraic expression to be devided
        denom (sympy Expr): Algebraic expression to divie.
            The Expr must be a single term.

    Returns:
        sympy Expr: result of division eq/denom.
    """
    return sp.Mul(eq, 1 / denom)
    # calculate Mul with evaluate=True,
    # because the result with evaluate=False will be like
    # `(1/(12*h))*(f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})`,
    # unlike the result we want as follows:
    # `(f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h)`
