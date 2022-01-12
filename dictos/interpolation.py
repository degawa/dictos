import sympy as sp
import utils as util


def getInterpolationEquation(stencil, sameSubscriptsAsStencil=False):
    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, util.DEFAULT_INTERVAL_SYMBOL_STR
    )
    fSet = util.create_set_of_function_symbols_at_coordinate(
        xSet, util.DEFAULT_FUNCTION_SYMBOL_STR, sameSubscriptsAsStencil
    )

    coef = getInterpolationCoefficients(stencil)

    return sp.simplify(sum([coef[i] * fSet[i] for i in range(len(fSet))]))


def getInterpolationCoefficients(stencil, as_numr_denom=False):
    import lagrangian_polynomial as lp

    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, util.DEFAULT_INTERVAL_SYMBOL_STR
    )
    fSet = util.create_set_of_function_symbols_at_coordinate(
        xSet, util.DEFAULT_FUNCTION_SYMBOL_STR
    )

    x = sp.symbols(util.DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR)

    doe = len(xSet) - 1
    eq = sum([lp.LagrangianBasis(x, doe, i, xSet) * fSet[i] for i in range(doe + 1)])

    num, den = sp.simplify(eq.subs(x, 0)).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [num / den_coef[0] for num in num_coef]

    return util.simplify_coefficients(coef, as_numr_denom)


def getTruncationError(stencil, intervalSymbolStr=util.DEFAULT_INTERVAL_SYMBOL_STR):
    import taylor_expansion as te

    xSet = util.create_set_of_coordinate_symbols_from_stencil(
        stencil, interval_symbol_str=intervalSymbolStr
    )

    coef = getInterpolationCoefficients(stencil)

    num_expterm = len(xSet)
    f_te = [te.TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i] * f_te[i] for i in range(len(xSet))])

    return sp.simplify(
        sp.symbols(util.DEFAULT_FUNCTION_SYMBOL_STR)
        - sp.nsimplify(eq, rational=True, tolerance=1e-10)
    ).as_leading_term()
