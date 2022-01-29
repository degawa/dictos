import sympy as sp
import numpy as np

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
