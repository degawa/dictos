import sympy as sp
import utils as util


def getFiniteDifferenceEquation(stencil, orderOfDifference=1,
                                intervalSymbolStr=util._DefaultIntervalSymbolStr,
                                sameSubscriptsAsStencil=False,
                                evaluate=True):
    xSet = util.createXSetFromStencil(stencil, intervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(xSet, util._DefaultFunctionSymbolStr,
                                                 sameSubscriptsAsStencil)

    if evaluate:
        coef = getFiniteDifferenceCoefficients(stencil, orderOfDifference)
        eq = sp.simplify(sum([coef[i]*fSet[i] for i in range(len(fSet))])
                         / sp.symbols(intervalSymbolStr)**orderOfDifference)
    else:
        numr, denom = getFiniteDifferenceCoefficients(
            stencil, orderOfDifference, as_numr_denom=True)
        eq = util.div(util.dotproduct_no_eval(numr, fSet), denom*sp.symbols(
            intervalSymbolStr)**orderOfDifference)

    return eq


def getFiniteDifferenceCoefficients(stencil, orderOfDifference=1, as_numr_denom=False):
    import lagrangianpoly as lp

    xSet = util.createXSetFromStencil(stencil, util._DefaultIntervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(
        xSet, util._DefaultFunctionSymbolStr)

    x = sp.symbols(util._DefaultIndependentVariableSymbolStr)
    num, den = lp.LagrangianPoly(x, xSet, fSet).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [lp.Derivative(num/den_coef[0], x, orderOfDifference)
            for num in num_coef]

    return util.simplifyCoefficients(coef, as_numr_denom)


def getTruncationError(stencil, orderOfDifference,
                       intervalSymbolStr=util._DefaultIntervalSymbolStr):
    import TaylorExpansion as te
    xSet = util.createXSetFromStencil(
        stencil, intervalSymbolStr=intervalSymbolStr)

    coef = getFiniteDifferenceCoefficients(xSet, orderOfDifference)

    num_expterm = len(xSet)+orderOfDifference
    f_te = [te.TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i]*f_te[i] for i in range(len(xSet))])

    intervalSymbol = sp.symbols(intervalSymbolStr)
    return sp.simplify(te._getDerivativeSymbol(util._DefaultFunctionSymbolStr,
                                               orderOfDifference)
                       - sp.nsimplify(eq/intervalSymbol**orderOfDifference,
                                      rational=True,
                                      tolerance=1e-10)).as_leading_term()
