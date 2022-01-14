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


def getFiniteDifferenceEquation(
    stencil,
    orderOfDifference=1,
    intervalSymbolStr=DEFAULT_INTERVAL,
    sameSubscriptsAsStencil=False,
    evaluate=True,
):
    xSet = create_coordinate_symbols(stencil, intervalSymbolStr)
    fSet = create_function_symbols(xSet, DEFAULT_FUNCTION, sameSubscriptsAsStencil)

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
            dot_product(numr, fSet),
            denom * sp.symbols(intervalSymbolStr) ** orderOfDifference,
        )

    return eq


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):
    xSet = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    fSet = create_function_symbols(xSet, DEFAULT_FUNCTION)

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)
    num, den = lagrangian_poly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [derivative(num / den_coef[0], x, orderOfDifference) for num in num_coef]

    return simplify_coefficients(coef, as_numr_denom)


def getTruncationError(stencil, orderOfDifference, intervalSymbolStr=DEFAULT_INTERVAL):
    xSet = create_coordinate_symbols(stencil, interval=intervalSymbolStr)

    coef = getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet) + orderOfDifference
    f_te = [taylor_series(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return sp.simplify(
        derivative_symbol(DEFAULT_FUNCTION, orderOfDifference)
        - sp.nsimplify(
            eq / intervalSymbol ** orderOfDifference, rational=True, tolerance=1e-10
        )
    ).as_leading_term()
