import sympy as sp
import numpy as np

from dictos.linalg import div
from dictos.stencil import get_subscript
from .error.internal import UnexpectedDenominatorError


def simplify_coefficients(coef, as_numer_denom: bool = False):
    """
    simplify coefficients in floating-point number
    to ratioanl numbers.
    When as_numer_denom flag is set, this returns the numberator and
    denominator separately.

    Args:
        coef (list of sympy Mul and numbers): finite difference and
            interpolation coefficients to be simplified.
        as_numer_denom (bool, optional): flag to return the numerator
            and denominator separately.
            Defaults to False.

    Returns:
        list of sympy Rational: simplified coefficients.
            or
        list of sympy numbers, int:
            numerator and denominator of coefficients.
            coefficients are commutative
            with the least common multiple of the denominator.
    """
    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]
    # extract numbers from list of coefficients like [1/h**2, ...].

    # rationalize coefficients and then divide those into numerator and denominator.
    # Execute while numbers other than one are contained in the coefficients' denominator.
    has_Rational = True
    max_denominator = 1000000  # sympy's default value
    significant_digits = 16  # approx significant digits of Float in decimal notation.
    while has_Rational and max_denominator <= 10 ** significant_digits:
        coef_rational = [
            sp.Rational(str(c)).limit_denominator(max_denominator) for c in coef_num
        ]
        # rationalize coefficients.
        # The argument of `sympy.Rational` must be a string in this case.
        # `limit_denominator` is used to avoid represent a real number
        # including the error coused by the binary representation.

        denom = [c.q for c in coef_rational]
        denom_lcm = int(np.lcm.reduce(np.array(denom)))
        # extract denomenator of each coefficient
        # and calculate the least common multiple.

        numer = [c * denom_lcm for c in coef_rational]
        # list of numerator divided by the least common multiple.

        has_Rational = any(n.q != 1 for n in numer)
        max_denominator *= 10
        # rationalization is incomplete
        # if the numerator contains at least a number other than 1.
        # retry rationalizaiton with larger maximum value of the denominator.

    if as_numer_denom:
        return numer, denom_lcm
    else:
        return [n / denom_lcm for n in numer]


def extract_coefficients_as_numer_denom(expr, f_set):
    """
    Extract coefficients as numerator and denominator
    from the polynomial in sympy Expr.

    Args:
        expr (sympy Expr): a polynomial from which coefficients are extracted.
        f_set (list or tuple of sympy symbols): set of functions.

    Raises:
        UnexpectedDenominatorError: if
            type of denom_coef is not list or tuple or
            type of denom_coef is list or tuple, but its length is greater than 1.

    Returns:
        list of sympy numbers, list [1] of sympy numbers:
            numerator and denominator of coefficients.
    """

    numer, denom = expr.as_numer_denom()
    numer_coef = numer.as_poly(f_set).coeffs()
    denom_coef = denom.as_poly(f_set).coeffs()
    # extract numerator and denominator from the polynomial

    def unexpected_type():
        return (type(denom_coef) is not list and type(denom_coef) is not tuple) or (
            (type(denom_coef) is list or type(denom_coef) is tuple)
            and len(denom_coef) > 1
        )

    if unexpected_type():
        raise UnexpectedDenominatorError(denom_coef)
        # raise error if
        # - type of denom_coef is not list or tuple
        # or
        # - type of denom_coef is list or tuple, but its length is greater than 1

    return numer_coef, denom_coef


def decompose_addition(expr_add):
    """
    Decompose additions into each term and extract all terms
    from sympy.core.add.Add expression.

    Args:
        expr_add (sympy.core.add.Add):
            an expression to extract all terms.
            for example, `0a + b + 2c` or `(0a + (b + 2c))`.

    Returns:
        list:
            all terms extracted from `expr_add`.
            The order of terms is not keeped
            when sympy.core.add.Add is contained in `expr_add`
            caused by a behavior of sympy.Add.make_args().

    Examples:
        >>> import sympy as sp
        >>> from sympy.abc import *
        >>> from dictos import utils as utl
        >>> numer = 1 * a + 2 * b + 3 * c
        >>> terms = utl.extract_terms(numer)
        >>> print(terms)
        [a, 2*b, 3*c]
        >>> numer = 0 * a + 1 * b + 2 * c
        >>> terms = utl.extract_terms(numer)
        >>> print(terms)
        [b, 2*c]
        >>> numer = sp.Mul(2, c, evaluate=False)
        >>> numer = sp.Add(numer, sp.Mul(1, b, evaluate=False), evaluate=False)
        >>> numer = sp.Add(numer, sp.Mul(0, a, evaluate=False), evaluate=False)
        >>> terms = utl.extract_terms(numer)
        >>> print(terms)
        [2*c, b, 0*a]
        >>> numer = sp.Mul(2, c, evaluate=False)
        >>> numer = sp.Add(sp.Mul(1, b, evaluate=False), numer, evaluate=False)
        >>> numer = sp.Add(sp.Mul(0, a, evaluate=False), numer, evaluate=False)
        >>> terms = extract_terms(numer)
        >>> print(terms)
        [0*a, b, 2*c]
    """

    if type(expr_add) is sp.core.symbol.Symbol or type(expr_add) is sp.core.mul.Mul:
        return expr_add
        # is type of `expr_add` is Symbol or Mul, there is nothing more to do.

    if type(expr_add) is not sp.core.add.Add:
        raise TypeError(expr_add, type(expr_add))
        # raise error if
        # - type of `expr_add` is not Add.

    terms = []
    # list for containing all terms.

    args = sp.Add.make_args(expr_add)
    # decompose and extract all terms from sympy.core.add.Add as tuple.
    # a + 2*b +3*c -> [a, 2*b, 3*c]

    if any(type(arg) is sp.core.add.Add for arg in args):
        # when sympy.core.add.Add is contained in args,
        # further decomposition is performed.
        for i in range(len(args)):
            if type(args[i]) is sp.core.add.Add:
                terms += decompose_addition(args[i])
                # when `args[i]` is sympy.core.add.Add,
                # all decomposed terms are added to the list.
            else:
                terms.append(args[i])
                # when `args[i]` is not sympy.core.add.Add,
                # add `args[i]` to the list as a term.
    else:
        terms = list(args[:])
        # when sympy.core.add.Add is not contained in args,
        # addition is successfully decomposed into each term.
        # return all terms as list.

    return terms


def sort_by_subscript(expr):
    """
    sort numerator of sympy expr by subscripts of symbols in the numerator.

    Args:
        expr (sympy expr):
            an expression to be sorted, like
            (-8*f_{-1} + f_{-2} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h),
            (4*f_{-1} - f_{-2} +  4*f_{1} - f_{2})/6

    Returns:
        sympy expr:
            a sorted expression, like
            (f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h),
            (-f_{-2} + 4*f_{-1} +  4*f_{1} - f_{2})/6
    """

    numer, denom = expr.as_numer_denom()
    # extract numerator and denominator from sympy expr.
    # `numer` is evaluated and a term multiplied by 0 is erased.
    # when expr is a fraction, `denom` is set to 1.

    # when denominator is an integer, values of `expr.as_numer_denom()`
    # and `expr.args` are different. for instance,
    # `((4*f_{-1} - f_{-2} +  4*f_{1} - f_{2})/6).as_numer_denom()`
    # returns `(4*f_{-1} - f_{-2} +  4*f_{1} - f_{2}, 6)`, but
    # `((4*f_{-1} - f_{-2} +  4*f_{1} - f_{2})/6).args` is
    # `(-f_{-2}/6, -f_{2}/6, 2*f_{-1}/3, 2*f_{1}/3)`.
    # so it is necessary to branch according to type of the denominator.

    if denom != 1 and type(denom) is not sp.core.numbers.Integer:
        # if `expr` is a fraction, and the denominator contains sympy symbol,
        numer = expr.args[-1]
        # extract the numerator.
        # `args` of assumed equation such as
        # (-8*f_{-1} + f_{-2} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h)
        # are
        # [1/12, 1/h, (-8*f_{-1} + f_{-2} + 0*f_{0} + 8*f_{1} - f_{2})].
        # so the numerator is expr.args[-1].

    terms = decompose_addition(numer)
    # decompose numerator into each term and extract all terms as a list.
    # intended result is [f_{-2}, -f_{2}, 0*f_{0}, -8*f_{-1}, 8*f_{1}].

    subscripts = []
    for n in terms:
        subscript = get_subscript(n)
        subscripts.append(float(subscript))
    # extract subscript from each term of the numerator.
    # intended result is [-2.000, 2.000, 0.000, -1.000, 1.000].

    numer_sorted_terms = [terms[i] for i in np.argsort(subscripts)]
    # sort terms of the numeartor by subscript.
    # a result of numpy.argsort(subscripts) is [0, 3, 2, 4, 1],
    # numer_sorted_terms is [terms[0], terms[3], terms[2], terms[4], terms[1]].

    # rearrange numerator in order of the subscripts.
    if all([s < 0 for s in subscripts]):
        # if all subscripts are negative integer,
        numer_sorted = numer_sorted_terms[0]
        for n in numer_sorted_terms[1:]:
            numer_sorted = sp.Add(numer_sorted, n, evaluate=False)
        # (((3*f_{-1}) - 3*f_{-2}) + f_{-3})
    elif all([s > 0 for s in subscripts]):
        # if all subscripts are positive integer,
        numer_sorted = numer_sorted_terms[-1]
        for n in numer_sorted_terms[-2::-1]:
            numer_sorted = sp.Add(numer_sorted, n, evaluate=False)
        # (((3*f_{1}) - 3*f_{2}) + f_{3})
    else:
        numer_sorted = numer_sorted_terms[-1]
        for n in numer_sorted_terms[-2::-1]:
            numer_sorted = sp.Add(numer_sorted, n, evaluate=False)
        # ((((-f_{2} + 8*f_{1}) + 0*f_{0}) - 8*f_{-1}) + f_{-2})

    if type(denom) is sp.core.numbers.Integer:
        return sp.Mul(numer_sorted, sp.Rational(1, denom), evaluate=False)
    else:
        return div(numer_sorted, denom)
    # reconstruct expr.
    # if the denominator contains a sympy symbol, use `div` function.
    # in the case that denominator is an integer, a result becomes
    # `-f_{-2}/6 + (4*f_{-1} + 4*f_{1} - f_{2})/6` when `div` is used.
    # to avoid this, an different expression is used.
