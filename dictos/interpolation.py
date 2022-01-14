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
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    fSet = create_function_symbols(xSet, DEFAULT_FUNCTION, same_subscripts_as_stencil)

    coef = coefficients(stencil)

    return sp.simplify(sum([coef[i] * fSet[i] for i in range(len(fSet))]))


def coefficients(stencil, as_numr_denom=False):
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    fSet = create_function_symbols(xSet, DEFAULT_FUNCTION)

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)

    doe = len(xSet) - 1
    eq = sum([lagrangian_basis(x, doe, i, xSet) * fSet[i] for i in range(doe + 1)])

    num, den = sp.simplify(eq.subs(x, 0)).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [num / den_coef[0] for num in num_coef]

    return simplify_coefficients(coef, as_numr_denom)


def truncation_error(stencil, interval=DEFAULT_INTERVAL):
    xSet = create_coordinate_symbols(stencil, interval=interval)

    coef = coefficients(stencil)

    num_expterm = len(xSet)
    f_te = [taylor_series(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    return sp.simplify(
        sp.symbols(DEFAULT_FUNCTION) - sp.nsimplify(eq, rational=True, tolerance=1e-10)
    ).as_leading_term()
