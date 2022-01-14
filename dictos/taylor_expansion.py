import sympy as sp

from .utils import DEFAULT_FUNCTION


def taylor_series(around, up_to):
    """
    calculate Taylor series of f(x) around h

    Args:
        around (sympy Mul): a position to calculate Taylor series of
            the function f.
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
    func = DEFAULT_FUNCTION
    # set the function symbol.
    # For a futere enhancement (specifying function symbol by users).

    deriv_orders = range(1, up_to + 1)
    # list of order of derivative [1, 2, 3, ..., up_to]

    df_set = [derivative_symbol(func, i) for i in deriv_orders]
    # derivatives of function f [f^(1), f^(2), f^(3), ..., f^(up_to)]

    h = around
    coef = [h ** i * sp.Rational(1, sp.factorial(i)) for i in deriv_orders]
    # coefficient each term [h, h**2/2, h**3/6, ..., h**up_to/up_to!]

    f = sp.symbols(func)
    series = f  # the first term
    for i in range(len(df_set)):
        series += df_set[i] * coef[i]
    # calculate summation of each term

    return series


def derivative_symbol(function, deriv):
    """
    returns n-th order derivative of function f as a sympy symbol

    Args:
        function (str): function symbol
        deriv (int): order of derivative

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
    return sp.symbols(function + "^(%d)" % deriv)
