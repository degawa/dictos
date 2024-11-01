import sympy as sp
import numpy as np

from dictos.linalg.linalg import div
from dictos.discrete.stencil import get_subscript
from dictos.utilities.exceptions.internal import UnexpectedDenominatorError


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

    Examples:
        >>> import sympy as sp
        >>> from sympy.abc import *
        >>> from dictos import utils as utl
        >>> expr = [1/(12*h), - 2/(3*h), + 2/(3*h), -1/(12*h)]
        >>> utl.simplify_coefficients(expr)
        [1/12, -2/3, 2/3, -1/12]
        >>> utl.simplify_coefficients(expr, as_numer_denom=True)
        ([1, -8, 8, -1], 12)
    """
    coef_num = [c if c.is_number else sp.poly(c).coeffs()[0] for c in coef]
    # extract numbers from list of coefficients like [1/h**2, ...].

    # convert each coefficient to a rational number
    coef_rational = [sp.nsimplify(c, rational=True) for c in coef_num]

    if not as_numer_denom:
        return coef_rational
    else:
        # Find the least common multiple of all denominators
        denom_lcm = sp.lcm([c.q for c in coef_rational])
        # extract denomenator of each coefficient
        # and calculate the least common multiple.

        # Multiply each coefficient by the LCM of denominators
        numer = [c * denom_lcm for c in coef_rational]

        return numer, denom_lcm


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
    from sympy.Add expression.

    Args:
        expr_add (sympy.Add):
            an expression to extract all terms.
            for example, `0a + b + 2c` or `(0a + (b + 2c))`.

    Returns:
        list:
            all terms extracted from `expr_add`.
            The order of terms is not keeped
            when sympy.Add is contained in `expr_add`
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

    if type(expr_add) is sp.Symbol or type(expr_add) is sp.Mul:
        return expr_add
        # if type of `expr_add` is Symbol or Mul, there is nothing more to do.

    if type(expr_add) is not sp.Add:
        raise TypeError(expr_add, type(expr_add))
        # raise error if
        # - type of `expr_add` is not Add.

    terms = []
    # list for containing all terms.

    args = sp.Add.make_args(expr_add)
    # decompose and extract all terms from sympy.Add as tuple.
    # a + 2*b +3*c -> [a, 2*b, 3*c]

    if any(type(arg) is sp.Add for arg in args):
        # when sympy.Add is contained in args,
        # further decomposition is performed.
        for i in range(len(args)):
            if type(args[i]) is sp.Add:
                terms += decompose_addition(args[i])
                # when `args[i]` is sympy.Add,
                # all decomposed terms are added to the list.
            else:
                terms.append(args[i])
                # when `args[i]` is not sympy.Add,
                # add `args[i]` to the list as a term.
    else:
        terms = list(args[:])
        # when sympy.Add is not contained in args,
        # addition is successfully decomposed into each term.
        # return all terms as list.

    return terms


def drop_coefficient_of_1(expr: sp.Expr) -> sp.Expr:
    """
    Drop the coefficient from terms with a coefficient of 1.

    Removes the coefficients of terms with a coefficient of 1 from the given expression.
    For example, `1*x + 2*y` is converted to `x + 2*y`.

    Args:
        expr (sp.Expr): sympy Expr. Symbol, Add, or Mul.

    Returns:
        sp.Expr: An expression from which the coefficient 1 is removed.
            The structure of the original expression is retained.

    Notes:
        - No evaluation of the original expression is performed (evaluate=False)
        - Coefficients of nested expressions, like `(1*x + y)*z`, are not processed correctly.
        - Depends on the `decompose_addition` method.

    Raises:
        TypeError:  if type of expr is not sympy Symbol, Add, or Mul.

    Examples:
        >>> from dictos.utilities import utils as utl
        >>> from sympy.abc import *
        >>> utl.drop_coefficient_of_1(1*x+2*y)
        x + 2*y
        >>> utl.drop_coefficient_of_1(2*x+3*y)
        2*x + 3*y
        >>> utl.drop_coefficient_of_1((1*x+2*y)*z)
        z
    """

    if type(expr) is sp.Symbol:
        return expr
        # if type of `expr` is Symbol, it is assumed to be unary with no coefficients
        # and returned as is.

    if type(expr) is sp.Add:
        args = decompose_addition(expr)
        # if expr is addition equation, decompose it into individual terms
    else:
        args = [expr]
        # for unary, convert to a list of single element for later summation using Add.

    def _drop_1(arg: sp.Mul) -> sp.Expr:
        """
        Remove a coefficient from a multiplication terms with a coefficeint of 1.

        Args:
            arg (sp.Mul): A multiplication term

        Returns:
            sp.Expr: The term with the coefficient removed if the coefficient is 1,
                otherwise the original term.
        """
        coef, terms = arg.as_coeff_mul()
        return terms[0] if coef == 1 else arg

    return sp.Add(
        *[_drop_1(arg) if type(arg) is sp.Mul else arg for arg in args], evaluate=False
    )
    # remove the coefficient of 1 from each term and reconstruct the addition equation.
    # each term is not evaluated (evaluate=False).


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

    if denom != 1 and type(denom) is not sp.Integer:
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
        #
        # 3*f_{-1} - 3*f_{-2} + f_{-3}
        numer_sorted = sp.Add(*numer_sorted_terms[::-1], evaluate=False)
    else:
        # 3*f_{1} - 3*f_{2} + f_{3}
        # f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2}
        # etc.
        numer_sorted = sp.Add(*numer_sorted_terms, evaluate=False)

    if type(denom) is sp.Integer:
        return sp.Mul(numer_sorted, sp.Rational(1, denom), evaluate=False)
    else:
        return div(numer_sorted, denom)
    # reconstruct expr.
    # if the denominator contains a sympy symbol, use `div` function.
    # in the case that denominator is an integer, a result becomes
    # `-f_{-2}/6 + (4*f_{-1} + 4*f_{1} - f_{2})/6` when `div` is used.
    # to avoid this, an different expression is used.
