import sympy as sp

import utils as util
import lagrangian_polynomial as lp
import taylor_expansion as te


def getFiniteDifferenceEquation(
    stencil,
    orderOfDifference=1,
    intervalSymbolStr=util.DEFAULT_INTERVAL_SYMBOL_STR,
    sameSubscriptsAsStencil=False,
    evaluate=True,
):
    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, intervalSymbolStr
    )
    fSet = util.create_set_of_function_symbols_at_coordinate(
        xSet, util.DEFAULT_FUNCTION_SYMBOL_STR, sameSubscriptsAsStencil
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
        eq = util.div(
            util.dotproduct(numr, fSet),
            denom * sp.symbols(intervalSymbolStr) ** orderOfDifference,
        )

    return eq


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):
    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, util.DEFAULT_INTERVAL_SYMBOL_STR
    )
    fSet = util.create_set_of_function_symbols_at_coordinate(
        xSet, util.DEFAULT_FUNCTION_SYMBOL_STR
    )

    x = sp.symbols(util.DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [lp.Derivative(num / den_coef[0], x, orderOfDifference) for num in num_coef]

    return util.simplify_coefficients(coef, as_numr_denom)


def getTruncationError(
    stencil, orderOfDifference, intervalSymbolStr=util.DEFAULT_INTERVAL_SYMBOL_STR
):
    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, interval_symbol_str=intervalSymbolStr
    )

    coef = getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet) + orderOfDifference
    f_te = [te.TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return sp.simplify(
        te._getDerivativeSymbol(util.DEFAULT_FUNCTION_SYMBOL_STR, orderOfDifference)
        - sp.nsimplify(
            eq / intervalSymbol ** orderOfDifference, rational=True, tolerance=1e-10
        )
    ).as_leading_term()
