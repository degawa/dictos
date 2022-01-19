import sympy as sp

from .utils import DEFAULT_FUNCTION
from .error.taylor_expansion import (
    UnsupportedOrderOfDerivativeError,
    NumberOfExpansionTermsIsNotNaturalNumberError,
)


def taylor_series(around, up_to):
    """
    calculate Taylor series of f(x) around h

    Args:
        around (sympy Mul): a relative position to calculate
            Taylor series of the function f.
            pass h not x+h when calculate Taylor series of f(x) at x+h.
        up_to (int): number of terms in the Taylor series
            excluding the first term

    Returns:
        sympy Expr: calculated Taylor series.

    Examples:
        >>> from dictos import taylor_expansion as te
        >>> import sympy as sp
        >>> h = sp.symbols("h")
        >>> te.taylor_series(h, 4)
        f + f^(1)*h + f^(2)*h**2/2 + f^(3)*h**3/6 + f^(4)*h**4/24
        >>> te.taylor_series(-h, 4)
        f - f^(1)*h + f^(2)*h**2/2 - f^(3)*h**3/6 + f^(4)*h**4/24
        >>> te.taylor_series(2*h, 4)
        f + 2*f^(1)*h + 2*f^(2)*h**2 + 4*f^(3)*h**3/3 + 2*f^(4)*h**4/3
    """
    num_term = up_to + 1
    # number of series expansion terms including the first term.
    if num_term < 0:
        raise NumberOfExpansionTermsIsNotNaturalNumberError(num_term)
        # raise error if
        # - number of series expansion terms is not the natural number.

    func = DEFAULT_FUNCTION
    # set the function symbol.
    # For a futere enhancement (specifying function symbol by users).

    deriv_orders = range(1, up_to + 1)  # +1 is correction for exclusive stop
    # list of order of derivative [1, 2, 3, ..., up_to]

    df_set = [derivative_symbol(func, i) for i in deriv_orders]
    # derivatives of function f [f^(1), f^(2), f^(3), ..., f^(up_to)]

    h = around
    coef = [h ** i * sp.Rational(1, sp.factorial(i)) for i in deriv_orders]
    # coefficient each term [h, h**2/2, h**3/6, ..., h**up_to/up_to!]

    f = sp.symbols(func)
    series = f + sum([df_set[i] * coef[i] for i in range(len(df_set))])
    # calculate summation of each term

    return series


def derivative_symbol(function, deriv):
    """
    returns n-th order derivative of function f as a sympy symbol

    Args:
        function (str): function symbol
        deriv (int): order of derivative

    Raises:
        UnsupportedOrderOfDerivativeError: if
            unsupproted order of derivative (< 0) is passed.

    Returns:
        sympy symbol: n-th derivative of function

    Examples:
        >>> from dictos import taylor_expansion as te
        >>> te.derivative_symbol("f", 1)
        f^(1)
        >>> te.derivative_symbol("f", 2)
        f^(2)
        >>> te.derivative_symbol("f", 3)
        f^(3)
    """
    if deriv < 0:
        raise UnsupportedOrderOfDerivativeError(deriv)
        # raise error if
        # - unsupproted order of derivative (< 0).

    return sp.symbols(function + "^(%d)" % deriv)
    # TODO: #33 return `"f"` not `"f^(0)"` when deriv is 0
