import sympy as sp

from .utils import DEFAULT_FUNCTION


def taylor_series(around, up_to):
    deriv_term = up_to
    df_set = sp.symbols(DEFAULT_FUNCTION + "^((1:{:d}))".format(deriv_term + 1))

    h = around
    coef = [h ** i * sp.Rational(1, sp.factorial(i)) for i in range(1, deriv_term + 1)]
    f = sp.symbols(DEFAULT_FUNCTION)
    series = f
    for i in range(len(df_set)):
        series += df_set[i] * coef[i]

    return series


def derivative_symbol(function, deriv):
    return sp.symbols(function + "^(%d)" % deriv)
