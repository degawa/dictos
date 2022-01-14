import sympy as sp

from .utils import (
    DEFAULT_INTERVAL,
    DEFAULT_FUNCTION,
    DEFAULT_INDEPENDENT_VARIABLE,
    create_coordinate_symbols,
    create_function_symbols,
    simplify_coefficients,
)
from .lagrangian_polynomial import lagrangian_basis
from .taylor_expansion import taylor_series


def equation(stencil, same_subscripts_as_stencil=False):
    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION, same_subscripts_as_stencil)

    coef = coefficients(stencil)

    return sp.simplify(sum([coef[i] * f_set[i] for i in range(len(f_set))]))


def coefficients(stencil, as_numr_denom=False):
    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION)

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)

    num_set = len(x_set)
    degree = num_set - 1
    eq = sum([lagrangian_basis(x, degree, i, x_set) * f_set[i] for i in range(num_set)])

    numr, denom = sp.simplify(eq.subs(x, 0)).as_numer_denom()
    numr_coef = numr.as_poly(f_set).coeffs()
    denom_coef = denom.as_poly(f_set).coeffs()

    coef = [num / denom_coef[0] for num in numr_coef]

    return simplify_coefficients(coef, as_numr_denom)


def truncation_error(stencil, interval=DEFAULT_INTERVAL):
    x_set = create_coordinate_symbols(stencil, interval=interval)

    coef = coefficients(stencil)

    num_term = len(x_set)
    f_te = [taylor_series(x, num_term) for x in x_set]

    eq = sum([coef[i] * f_te[i] for i in range(len(x_set))])

    return sp.expand(
        sp.simplify(
            sp.symbols(DEFAULT_FUNCTION)
            - sp.nsimplify(eq, rational=True, tolerance=1e-10)
        )
    ).as_leading_term()
