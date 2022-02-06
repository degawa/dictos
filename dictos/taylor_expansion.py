import sympy as sp

from .spec import DEFAULT_DIFFERENTIAND, is_not_natrual_number
from .linalg import dot_product
from .error.taylor_expansion import (
    UnsupportedOrderOfDerivativeError,
    NumberOfExpansionTermsIsNotNaturalNumberError,
)


def taylor_series(around, up_to: int):
    """
    calculate Taylor series of f(x) around h

    Args:
        around (sympy Mul): a relative position to calculate
            Taylor series of the function f.
            pass h not x+h when calculate Taylor series of f(x) at x+h.
        up_to (int): number of terms in the Taylor series
            excluding the first term

    Raises:
        NumberOfExpansionTermsIsNotNaturalNumberError: if
            number of series expansion terms including the first term
            is not the natural number.

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
    if is_not_natrual_number(num_term):
        raise NumberOfExpansionTermsIsNotNaturalNumberError(num_term)
        # raise error if
        # - number of series expansion terms is not the natural number.

    func = DEFAULT_DIFFERENTIAND
    # set the function symbol.
    # For a futere enhancement (specifying function symbol by users).

    deriv_orders = range(up_to + 1)  # +1 is correction for exclusive stop
    # list of order of derivative [0, 1, 2, 3, ..., up_to]

    df_set = [derivative_symbol(func, i) for i in deriv_orders]
    # derivatives of function f [f, f^(1), f^(2), f^(3), ..., f^(up_to)]

    h = around
    coef = [h ** i * sp.Rational(1, sp.factorial(i)) for i in deriv_orders]
    # coefficient each term [1, h, h**2/2, h**3/6, ..., h**up_to/up_to!]

    series = dot_product(df_set, coef)
    # calculate summation of each term

    return series


def derivative_symbol(function: str, deriv: int):
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
    if is_not_natrual_number(deriv, include_zero=True):
        raise UnsupportedOrderOfDerivativeError(deriv)
        # raise error if
        # - unsupproted order of derivative (< 0).

    return sp.symbols(function) if deriv == 0 else sp.symbols(function + f"^({deriv})")
    # when deriv == 0, return `"f"` not `"f^(0)"`
    # when deriv >  0, return `"f^(deriv)"`
