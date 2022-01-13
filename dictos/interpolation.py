import sympy as sp

from .utils import (
    DEFAULT_INTERVAL_SYMBOL_STR,
    DEFAULT_FUNCTION_SYMBOL_STR,
    DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR,
    create_coordinate_symbols,
    create_set_of_function_symbols_at_coordinate,
    simplify_coefficients,
)
from .lagrangian_polynomial import LagrangianBasis
from .taylor_expansion import TaylorExpansion


def getInterpolationEquation(stencil, sameSubscriptsAsStencil=False):
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL_SYMBOL_STR)
    fSet = create_set_of_function_symbols_at_coordinate(
        xSet, DEFAULT_FUNCTION_SYMBOL_STR, sameSubscriptsAsStencil
    )

    coef = getInterpolationCoefficients(stencil)

    return sp.simplify(sum([coef[i] * fSet[i] for i in range(len(fSet))]))


def getInterpolationCoefficients(stencil, as_numr_denom=False):
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL_SYMBOL_STR)
    fSet = create_set_of_function_symbols_at_coordinate(
        xSet, DEFAULT_FUNCTION_SYMBOL_STR
    )

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR)

    doe = len(xSet) - 1
    eq = sum([LagrangianBasis(x, doe, i, xSet) * fSet[i] for i in range(doe + 1)])

    num, den = sp.simplify(eq.subs(x, 0)).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [num / den_coef[0] for num in num_coef]

    return simplify_coefficients(coef, as_numr_denom)


def getTruncationError(stencil, intervalSymbolStr=DEFAULT_INTERVAL_SYMBOL_STR):
    xSet = create_coordinate_symbols(stencil, interval=intervalSymbolStr)

    coef = getInterpolationCoefficients(stencil)

    num_expterm = len(xSet)
    f_te = [TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    return sp.simplify(
        sp.symbols(DEFAULT_FUNCTION_SYMBOL_STR)
        - sp.nsimplify(eq, rational=True, tolerance=1e-10)
    ).as_leading_term()
