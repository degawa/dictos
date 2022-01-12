import sympy as sp
import numpy as np
import fractions as fr


DEFAULT_INDEPENDENT_VARIABLE_SYMBOL_STR = "x"
DEFAULT_INTERVAL_SYMBOL_STR = "h"
DEFAULT_FUNCTION_SYMBOL_STR = "f"


def create_set_of_coordinate_symbols_from_stencil(
    stencil, interval_symbol_str=DEFAULT_INTERVAL_SYMBOL_STR
):
    # todo raise error when len(stencil)==0
    return [stencil[i] * sp.symbols(interval_symbol_str) for i in range(len(stencil))]


def create_set_of_function_symbols_at_coordinate(
    x_set,
    function_symbol_str=DEFAULT_FUNCTION_SYMBOL_STR,
    same_subscripts_as_stencil=False,
):

    if same_subscripts_as_stencil:
        stencil = [x if x.is_number else sp.poly(x).coeffs()[0] for x in x_set]
        subscript = ["_{%d}" % i if i.is_integer else "_{%2.1f}" % i for i in stencil]
        str = "".join([function_symbol_str + s + " " for s in subscript])
        fSet = sp.symbols(str)
    else:
        fSet = sp.symbols(function_symbol_str + "_0:{:d}".format(len(x_set)))

    return fSet


def simplify_coefficients(coef, as_numr_denom):
    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]

    coef_rational = [
        sp.Rational(fr.Fraction(str(c)).limit_denominator(100000)) for c in coef_num
    ]

    denom = [c.q for c in coef_rational]
    denom_lcm = np.lcm.reduce(np.array(denom))
    numr = [c * denom_lcm for c in coef_rational]

    if as_numr_denom:
        return numr, denom_lcm
    else:
        return [n / denom_lcm for n in numr]


def dotproduct(numr, f_set):
    idx_start = len(f_set) - 1
    idx_end = -1
    eq = sp.Mul(numr[idx_start], f_set[idx_start], evaluate=False)
    for i in range(idx_start - 1, idx_end, -1):
        eq = sp.Add(eq, sp.Mul(numr[i], f_set[i], evaluate=False), evaluate=False)

    return eq


def div(eq, denom):
    return sp.Mul(eq, 1 / denom)
