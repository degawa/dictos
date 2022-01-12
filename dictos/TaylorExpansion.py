import sympy as sp
import utils as util


def TaylorExpansion(h, n):
    df_set = sp.symbols(util._DefaultFunctionSymbolStr +
                        '^((1:{:d}))'.format(n+1))

    coef = [h**i*sp.Rational(1, sp.factorial(i)) for i in range(1, n+1)]
    f = sp.symbols(util._DefaultFunctionSymbolStr)
    te = f
    for i in range(len(df_set)):
        te += df_set[i]*coef[i]

    return te


def _getDerivativeSymbol(functionSymbolStr, n):
    return sp.symbols(functionSymbolStr+"^(%d)" % n)
