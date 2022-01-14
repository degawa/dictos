import sympy as sp

from .utils import DEFAULT_FUNCTION


def taylor_series(around, up_to):
    func = DEFAULT_FUNCTION

    deriv_orders = range(1, up_to + 1)
    df_set = [derivative_symbol(func, i) for i in deriv_orders]

    h = around
    coef = [h ** i * sp.Rational(1, sp.factorial(i)) for i in deriv_orders]

    f = sp.symbols(func)
    series = f
    for i in range(len(df_set)):
        series += df_set[i] * coef[i]

    return series


def derivative_symbol(function, deriv):
    return sp.symbols(function + "^(%d)" % deriv)
