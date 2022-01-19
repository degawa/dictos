import sympy as sp

from .utils import (
    DEFAULT_INTERVAL,
    DEFAULT_FUNCTION,
    DEFAULT_INDEPENDENT_VARIABLE,
    create_coordinate_symbols,
    create_function_symbols,
    simplify_coefficients,
    extract_coefficients_as_numer_denom,
    has_zero,
)
from .lagrangian_polynomial import lagrangian_basis, lagrangian_poly
from .taylor_expansion import taylor_series
from .error.stencil import ContainsZeroError


def equation(stencil, same_subscripts_as_stencil=False):
    """
    derive interpolation equation based on given stencil.
    The equation compute interpolation at stencil = 0

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        same_subscripts_as_stencil (bool, optional): flag
            to make function subscripts the same as the stencil.
            Defaults to False, the subscripts start from 0.

    Returns:
        sympy Expr: derived interpolation equation.

    Examples:
        >>> from dictos import interpolation as intp
        >>> intp.equation([-1, 1])
        f_0/2 + f_1/2
        >>> intp.equation([-1.5, -0.5, 0.5, 1.5])
        -f_0/16 + 9*f_1/16 + 9*f_2/16 - f_3/16
        >>> intp.equation([-1.5, -0.5, 0.5, 1.5],same_subscripts_as_stencil=True)
        9*f_{-0.5}/16 - f_{-1.5}/16 + 9*f_{0.5}/16 - f_{1.5}/16
    """

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION, same_subscripts_as_stencil)
    # create set of coordinate and function symbols from stencil.
    # [-2, -1, 1, 2] -> [-2*h, -h, h, 2*h]
    # [-2*h, -h, h, 2*h] -> [f_0, f_1, f_2, f_3]
    #                       [f_{-2}, f_{-1}, f_{1}, f_{2}]
    #                       (same_subscripts_as_stencil = True)

    coef = coefficients(stencil)
    # derive interpolation coefficients based on given stencil

    return sp.simplify(sum([coef[i] * f_set[i] for i in range(len(f_set))]))
    # calculate dot product of the coef and function values.


def coefficients(stencil, as_numr_denom=False):
    """
    derive interpolation coefficients based on given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        as_numr_denom (bool, optional): flag to return the numerator
            and denominator separately.
            Defaults to False.

    Raises:
        ContainsZeroError: if stencil contains 0.

    Returns:
        list of sympy Rational: simplified coefficients.
            or
        list of sympy numbers, int:
            numerator and denominator of coefficients.
            coefficients are commutative
            with the least common multiple of the denominator.

    Examples:
        >>> from dictos import interpolation as intp
        >>> intp.coefficients([-1, 1])
        [1/2, 1/2]
        >>> intp.coefficients([-1.5, -0.5, 0.5, 1.5])
        [-1/16, 9/16, 9/16, -1/16]
        >>> intp.coefficients([-1.5, -0.5, 0.5, 1.5], as_numr_denom=True)
        ([-1, 9, 9, -1], 16)
    """
    if has_zero(stencil):
        raise ContainsZeroError
        # raise error if stencil contains 0

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION)
    # create set of coordinate and function symbols from stencil.
    # [-2, -1, 1, 2] -> [-2*h, -h, h, 2*h]
    # [-2*h, -h, h, 2*h] -> [f_0, f_1, f_2, f_3]

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)
    eq = lagrangian_poly(x, x_set, f_set)
    # construct lagrangian interpolation polynomial

    poly = sp.simplify(eq.subs(x, 0))
    numr_coef, denom_coef = extract_coefficients_as_numer_denom(poly, f_set)
    # extract numerator and denomitaor the polynomial

    coef = [num / denom_coef[0] for num in numr_coef]
    # get coefficients of each terms as a list. Another expression
    # `derivative(lagrangian_poly(x, x_set, f_set), x, deriv).as_poly(f_set).coeffs()`
    # erases terms with a coefficient of 0.

    return simplify_coefficients(coef, as_numr_denom)
    # simplify floating-point number coefficients to ratioanl numbers


def truncation_error(stencil, interval=DEFAULT_INTERVAL):
    """
    derive the leading-order of error term
    in the interpolation equation based on the given stencil

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL.

    Returns:
        sympy Expr: the leading-order of error term

    Examples
        >>> from dictos import interpolation as intp
        >>> intp.truncation_error([-1, 1])
        -f^(2)*h**2/2
        >>> intp.truncation_error([-1.5, -0.5, 0.5, 1.5])
        3*f^(4)*h**4/128
    """

    coef = coefficients(stencil)
    # derive interpolation coefficients based on given stencil

    x_set = create_coordinate_symbols(stencil, interval=interval)
    # create set of coordinate symbols from stencil.
    # [-2, -1, 1, 2] -> [-2*h, -h, h, 2*h]

    num_term = len(x_set)
    f_te = [taylor_series(x, num_term) for x in x_set]
    # calculate Taylor series around points in x_set.

    intp_eq = sum([coef[i] * f_te[i] for i in range(len(x_set))])
    # calculate weighted sum of Taylor series.
    # for instance, 2nd-order 2-point interpolation is
    # intp_eq [= f(h)/2 + f(-h)/2)] = f(0) + f^(2)*h**2/2 + ...

    h = sp.symbols(interval)
    return sp.expand(
        sp.simplify(
            sp.symbols(DEFAULT_FUNCTION)
            - sp.nsimplify(intp_eq, rational=True, tolerance=1e-10)
        )
    ).as_leading_term(h)
    # extract the leading-order of errer term.
    # A interpolation formulation with error term is, for instance,
    # f(0) = (f(h) + f(-h))/2 - f^(2)*h**2/2 - ...
    # to extract error terms, reformulate intp_eq as
    # f(0) - intp_eq = - f^(2)*h**2/2 - ...
