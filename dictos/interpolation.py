import sympy as sp
import utils as util


def getInterpolationEquation(stencil,
                             sameSubscriptsAsStencil=False):
    xSet = util.createXSetFromStencil(stencil, util._DefaultIntervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(xSet, util._DefaultFunctionSymbolStr,
                                                 sameSubscriptsAsStencil)

    coef = getInterpolationCoefficients(stencil)

    return sp.simplify(sum([coef[i]*fSet[i] for i in range(len(fSet))]))


def getInterpolationCoefficients(stencil, as_numr_denom=False):
    import lagrangianpoly as lp

    xSet = util.createXSetFromStencil(stencil, util._DefaultIntervalSymbolStr)
    fSet = util.createSetOfFunctionSymbolsAtXSet(
        xSet, util._DefaultFunctionSymbolStr)

    x = sp.symbols(util._DefaultIndependentVariableSymbolStr)

    doe = len(xSet)-1
    eq = sum([lp.LagrangianBasis(x, doe, i, xSet)*fSet[i]
              for i in range(doe+1)])

    num, den = sp.simplify(eq.subs(x, 0)).as_numer_denom()
    num_coef = num.as_poly(fSet).coeffs()
    den_coef = den.as_poly(fSet).coeffs()

    coef = [num/den_coef[0] for num in num_coef]

    return util.simplifyCoefficients(coef, as_numr_denom)


def getTruncationError(stencil,
                       intervalSymbolStr=util._DefaultIntervalSymbolStr):
    import TaylorExpansion as te
    xSet = util.createXSetFromStencil(
        stencil, intervalSymbolStr=intervalSymbolStr)

    coef = getInterpolationCoefficients(stencil)

    num_expterm = len(xSet)
    f_te = [te.TaylorExpansion(x, num_expterm) for x in xSet]

    eq = sum([coef[i]*f_te[i] for i in range(len(xSet))])

    return sp.simplify(sp.symbols(util._DefaultFunctionSymbolStr)
                       - sp.nsimplify(eq,
                                      rational=True,
                                      tolerance=1e-10)).as_leading_term()
