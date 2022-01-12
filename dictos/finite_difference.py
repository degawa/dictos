import sympy as sp

from .utils import (
    DEFAULT_FUNCTION_SYMBOL_STR,
    DEFAULT_INTERVAL_SYMBOL_STR,
    DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR,
    create_coordinate_symbols,
    create_set_of_function_symbols_at_coordinate,
    simplify_coefficients,
    dotproduct,
    div,
)
from .lagrangian_polynomial import LagrangianPoly, Derivative
from .taylor_expansion import TaylorExpansion, _getDerivativeSymbol


def getFiniteDifferenceEquation(
    stencil,
    orderOfDifference=1,
    intervalSymbolStr=DEFAULT_INTERVAL_SYMBOL_STR,
    sameSubscriptsAsStencil=False,
    evaluate=True,
):
    xSet = create_coordinate_symbols(stencil, intervalSymbolStr)
    fSet = create_set_of_function_symbols_at_coordinate(
        xSet, DEFAULT_FUNCTION_SYMBOL_STR, sameSubscriptsAsStencil
    )

    if evaluate:
        coef = getFiniteDifferenceCoefficients(stencil, orderOfDifference)
        eq = sp.simplify(
            sum([coef[i] * fSet[i] for i in range(len(fSet))])
            / sp.symbols(intervalSymbolStr) ** orderOfDifference
        )
    else:
        numr, denom = getFiniteDifferenceCoefficients(
            stencil, orderOfDifference, as_numr_denom=True
        )
        eq = div(
            dotproduct(numr, fSet),
            denom * sp.symbols(intervalSymbolStr) ** orderOfDifference,
        )

    return eq


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL_SYMBOL_STR)
    fSet = create_set_of_function_symbols_at_coordinate(
        xSet, DEFAULT_FUNCTION_SYMBOL_STR
    )

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR)
    num, den = LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [Derivative(num / den_coef[0], x, orderOfDifference) for num in num_coef]

    return simplify_coefficients(coef, as_numr_denom)


def getTruncationError(
    stencil, orderOfDifference, intervalSymbolStr=DEFAULT_INTERVAL_SYMBOL_STR
):
    xSet = create_coordinate_symbols(stencil, interval_symbol_str=intervalSymbolStr)

    coef = getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet) + orderOfDifference
    f_te = [TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return sp.simplify(
        _getDerivativeSymbol(DEFAULT_FUNCTION_SYMBOL_STR, orderOfDifference)
        - sp.nsimplify(
            eq / intervalSymbol ** orderOfDifference, rational=True, tolerance=1e-10
        )
    ).as_leading_term()
