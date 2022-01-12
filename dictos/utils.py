import sympy as sp

_DefaultIndependentVariableSymbolStr = 'x'
_DefaultIntervalSymbolStr = 'h'
_DefaultFunctionSymbolStr = 'f'


def createXSetFromStencil(stencil, intervalSymbolStr=_DefaultIntervalSymbolStr):
    return [stencil[i]*sp.symbols(intervalSymbolStr) for i in range(len(stencil))]


def createSetOfFunctionSymbolsAtXSet(xSet, functionSymbolStr=_DefaultFunctionSymbolStr,
                                     sameSubscriptsAsStencil=False):

    if sameSubscriptsAsStencil:
        stencil = [x if x.is_number else sp.poly(x).coeffs()[0] for x in xSet]
        subscript = ['_{%d}' % i if i.is_integer else '_{%2.1f}' % i
                     for i in stencil]
        str = ''.join([functionSymbolStr + s + ' ' for s in subscript])
        fSet = sp.symbols(str)
    else:
        fSet = sp.symbols(functionSymbolStr+'_0:{:d}'.format(len(xSet)))

    return fSet


def simplifyCoefficients(coef, as_numr_denom):
    import numpy as np
    import fractions as fr

    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]

    coef_rational = [sp.Rational(fr.Fraction(
        str(c)).limit_denominator(100000)) for c in coef_num]

    denom = [c.q for c in coef_rational]
    denom_lcm = np.lcm.reduce(np.array(denom))
    numr = [c*denom_lcm for c in coef_rational]

    if as_numr_denom:
        return numr, denom_lcm
    else:
        return [n/denom_lcm for n in numr]


def dotproduct_no_eval(numr, fSet):
    idx_start = len(fSet)-1
    idx_end = -1
    eq = sp.Mul(numr[idx_start], fSet[idx_start], evaluate=False)
    for i in range(idx_start-1, idx_end, -1):
        eq = sp.Add(eq, sp.Mul(numr[i], fSet[i],
                               evaluate=False), evaluate=False)

    return eq


def div(eq, denom):
    return sp.Mul(eq, 1/denom)
