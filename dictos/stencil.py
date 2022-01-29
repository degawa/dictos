import sympy as sp

from .spec import (
    DEFAULT_INTERVAL,
    DEFAULT_FUNCTION,
    has_duplicated_points,
    narrower_than_minimum_width,
)
from .error.stencil import TooNarrowError, DuplicatedPointError


def create_coordinate_symbols(stencil: list, interval: str = DEFAULT_INTERVAL) -> list:
    """
    create set of coordinate symbols from stencil.
    input a list of numbers like `[-1, 0, 1]` as stencil,
    this returns a list of coordinates like `[-h, 0, h]`.

    Args:
        stencil (list of int or float): stencil on regular or
            staggered grid. It is not allowed that a number
            in the list appears more than once.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL.

    Raises:
        TooNarrowError: if stencil is too narrow.
        DuplicatedPointError: if at least a number in the stencil
            appears more than once.

    Returns:
        list of sympy symbols: list of coordinates
            corresponding to the stencil.
    """

    if narrower_than_minimum_width(stencil):
        raise TooNarrowError(stencil)
        # raise error if
        # - stencil is too narrow to coompute finite difference or interpolation
    if has_duplicated_points(stencil):
        raise DuplicatedPointError(stencil)
        # raise error if
        # - at least a number in the stencil appears more than once.

    sorted_stencil = stencil[:]
    sorted_stencil.sort()
    # sorting the stencil to obtain an equation
    # in which terms are arranged in the order of the stencil.
    # a local variable is used to not change the passed stencil.

    return [sorted_stencil[i] * sp.symbols(interval) for i in range(len(stencil))]


def create_function_symbols(
    x_set: list,
    function: str = DEFAULT_FUNCTION,
    same_subscripts_as_stencil: bool = False,
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

        subscript = [f"{i}" if i.is_integer else f"{float(i):2.1f}" for i in stencil]
        str = "".join([function + "_{" + s + "}" + " " for s in subscript])
        # make string "f_{-1} f_{-0.5} ..." to pass `sympy.symbols`.
        # tail space is ignored in `sympy.symbols`

        f_set = sp.symbols(str)
        if type(f_set) == sp.core.symbol.Symbol:
            f_set = (f_set,)
        # make the return value's type a tuple.
        # sp.symbols(str) returns sp.core.symbol.Symbol
        # when str does not contain space
        # althought sp.symbols(f"_0:{n}") returns tuple
        # when n==1.
    else:
        f_set = sp.symbols(function + f"_0:{len(x_set)}")
        # make a tuple of sympy symbols from string.

    return f_set