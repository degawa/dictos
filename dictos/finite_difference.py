import sympy as sp

from .utils import (
    DEFAULT_FUNCTION,
    DEFAULT_INTERVAL,
    DEFAULT_INDEPENDENT_VARIABLE,
    create_coordinate_symbols,
    create_function_symbols,
    simplify_coefficients,
    dot_product,
    div,
)
from .lagrangian_polynomial import lagrangian_poly, derivative
from .taylor_expansion import taylor_series, derivative_symbol


def equation(
    stencil,
    deriv=1,
    interval=DEFAULT_INTERVAL,
    same_subscripts_as_stencil=False,
    evaluate=True,
):
    x_set = create_coordinate_symbols(stencil, interval)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION, same_subscripts_as_stencil)

    if evaluate:
        coef = coefficients(stencil, deriv)
        eq = sp.simplify(
            sum([coef[i] * f_set[i] for i in range(len(f_set))])
            / sp.symbols(interval) ** deriv
        )
    else:
        numr, denom = coefficients(stencil, deriv, as_numr_denom=True)
        eq = div(
            dot_product(numr, f_set),
            denom * sp.symbols(interval) ** deriv,
        )

    return eq


def coefficients(stencil, deriv=1, as_numr_denom=False):
    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION)

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)
    numr, denom = lagrangian_poly(x, x_set, f_set).as_numer_denom()
    numr_coef = numr.as_poly(f_set).coeffs()
    denom_coef = denom.as_poly(f_set).coeffs()

    coef = [derivative(num / denom_coef[0], x, deriv) for num in numr_coef]

    return simplify_coefficients(coef, as_numr_denom)


def truncation_error(stencil, deriv, interval=DEFAULT_INTERVAL):
    x_set = create_coordinate_symbols(stencil, interval=interval)

    coef = coefficients(x_set, deriv)

    num_term = len(x_set) + deriv
    f_ts = [taylor_series(x, num_term) for x in x_set]

    eq = sum([coef[i] * f_ts[i] for i in range(len(x_set))])

    h = sp.symbols(interval)
    return sp.simplify(
        derivative_symbol(DEFAULT_FUNCTION, deriv)
        - sp.nsimplify(eq / h ** deriv, rational=True, tolerance=1e-10)
    ).as_leading_term()
