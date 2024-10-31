import sympy as sp

from dictos.utilities.spec import (
    DEFAULT_INTERVAL,
    DEFAULT_DIFFERENTIAND,
    has_duplicated_points,
    narrower_than_minimum_width,
)
from dictos.error.stencil import TooNarrowError, DuplicatedPointError


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

    return [s * sp.symbols(interval) for s in sorted_stencil]


def create_differentiand_symbols(
    x_set: list,
    differentiand: str = DEFAULT_DIFFERENTIAND,
):
    """
    create set of differentiand symbols at coordinates from x_set.
    input a list of sympy symbols like `[-h, 0, h]` as x_set,
    this returns a list of differentiands at x_set like
    `[f_{-1}, f_{0}, f_{1}]`.

    Args:
        x_set (list of sympy symbols): coordinates on regular or
            staggered grid corresponding to the stencil.
        differentiand (str, optional): a differentiand symbol like `f`.
            Defaults to DEFAULT_DIFFERENTIAND.

    Returns:
        tuple of sympy symbols: tuple of differentiands at passed coordinates.
            Subscripts are enclused in curly braces like `f_{0}`.
    """
    # make differentiand subscripts the same as the stencil

    stencil = [x if x.is_number else sp.poly(x).coeffs()[0] for x in x_set]
    # extract numbers from list of coordinates.
    # coordinate consists of a number and a symbol,
    # such as -2*h and 1.5*h, extract the number as coefficint.
    # coordinate is 0, that is a number, use 0 as the stencil

    subscript = [to_subscript(i) for i in stencil]
    str = "".join([differentiand + "_{" + s + "}" + " " for s in subscript])
    # make string "f_{-1} f_{-0.5} ..." to pass `sympy.symbols`.
    # tail space is ignored in `sympy.symbols`

    f_set = sp.symbols(str)
    if type(f_set) == sp.Symbol:
        f_set = (f_set,)
    # make the return value's type a tuple.
    # sp.symbols(str) returns sp.Symbol
    # when str does not contain space
    # althought sp.symbols(f"_0:{n}") returns tuple
    # when n==1.

    return f_set


def to_subscript(number):
    """
    convert int or float to stencil subscript.

    Args:
        number (int, float, or sympy Number): a number to be a subscript.

    Returns:
        str: converted subscript.
    """
    n = sp.Rational(number)
    return f"{n}" if n.is_integer else f"{float(n):2.1f}"
    # if n is a float number, converted to "xx.x"
    # 2-digit means the maximum stencil width is 99.
    # 1 digimal place is enough because equidistance grid is supported.


def get_subscript(a_term):
    """
    get subscript from a differentiand symbol.

    Args:
        a_term (sympy Symbol or Mul):
            a differentiand symbol with subscript like `27*f_{-1}`

    Raises:
        TypeError: if `a_term` is not sympy Symbol or Mul.

    Returns:
        str: string of subscript in a differentiand symbol, like "-1"
    """

    if type(a_term) is sp.Symbol:
        # if `a_term` is a symbol with subscript and without coefficient
        # like `f_{-1}`, `f = f_{-1}`.
        f = a_term
    elif type(a_term) is sp.Mul:
        # if `a_term` is a symbol with subscript and coefficient
        # like `27*f_{-1}`, `f = f_{-1}`.
        i = 1 if type(a_term.args[1]) is sp.Symbol else 0
        f = a_term.args[i]
        # if `a_term.args` is (27, f_{-1}), chose `f_{-1}`.
    else:
        raise TypeError(a_term, type(a_term))
        # raise error if
        # - `a_term` is not sympy Symbol or Mul.

    f_str = str(f)
    start = f_str.find("{") + 1
    stop = f_str.find("}")
    subscript = f_str[start:stop]
    # convert a symbol `f_{-1}` to a string "f_{-1}",
    # then find potisions at "{" and "}", finally extract subscript
    # between "{" and "}".

    return subscript
