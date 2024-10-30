import sympy as sp

from .spec import (
    DEFAULT_INTERVAL,
    DEFAULT_DIFFERENTIAND,
    DEFAULT_INDEPENDENT_VARIABLE,
    has_zero,
)
from .stencil import create_coordinate_symbols, create_differentiand_symbols
from .utils import (
    simplify_coefficients,
    extract_coefficients_as_numer_denom,
    sort_by_subscript,
)
from .linalg import dot_product
from .lagrangian_polynomial import lagrangian_poly
from .taylor_expansion import taylor_series
from .error.stencil import ContainsZeroError
from .core.expr import Expr


def equation(stencil: list, sort: bool = True) -> Expr:
    """
    derive interpolation equation based on given stencil.
    The equation compute interpolation at stencil = 0

    Args:
        stencil (list of int): relative point numbers
            used for discretization.

    Returns:
        sympy Expr: derived interpolation equation.

    Examples:
        >>> from dictos import interpolation as intp
        >>> intp.equation([-1, 1])
        f_{-1}/2 + f_{1}/2
        >>> intp.equation([-1.5, -0.5, 0.5, 1.5])
        (-f_{-1.5} + 9*f_{-0.5} + 9*f_{0.5} - f_{1.5})/16
    """

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_differentiand_symbols(x_set, DEFAULT_DIFFERENTIAND)
    # create set of coordinate and differentiand symbols from stencil.
    # [-2, -1, 1, 2] -> [-2*h, -h, h, 2*h]
    # [-2*h, -h, h, 2*h] -> [f_{-2}, f_{-1}, f_{1}, f_{2}]

    coef = coefficients(stencil)
    # derive interpolation coefficients based on given stencil

    eq = dot_product(coef, f_set)
    # calculate dot product of the coefs and differentiands.

    if sort:
        eq = sort_by_subscript(eq)
    # sort numerator by subscript.

    return Expr(eq)


def coefficients(stencil: list, as_numer_denom: bool = False):
    """
    derive interpolation coefficients based on given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        as_numer_denom (bool, optional): flag to return the numerator
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
        >>> intp.coefficients([-1.5, -0.5, 0.5, 1.5], as_numer_denom=True)
        ([-1, 9, 9, -1], 16)
    """
    if has_zero(stencil):
        raise ContainsZeroError
        # raise error if stencil contains 0

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_differentiand_symbols(x_set, DEFAULT_DIFFERENTIAND)
    # create set of coordinate and differentiand symbols from stencil.
    # [-2, -1, 1, 2] -> [-2*h, -h, h, 2*h]
    # [-2*h, -h, h, 2*h] -> [f_{-2}, f_{-1}, f_{1}, f_{2}]

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)
    eq = lagrangian_poly(x, x_set, f_set)
    # construct lagrangian interpolation polynomial

    poly = sp.simplify(eq.subs(x, 0))
    numer_coef, denom_coef = extract_coefficients_as_numer_denom(poly, f_set)
    # extract numerator and denomitaor the polynomial

    coef = [num / denom_coef[0] for num in numer_coef]
    # get coefficients of each terms as a list. Another expression
    # `derivative(lagrangian_poly(x, x_set, f_set), x, deriv).as_poly(f_set).coeffs()`
    # erases terms with a coefficient of 0.

    return simplify_coefficients(coef, as_numer_denom)
    # simplify floating-point number coefficients to ratioanl numbers


def truncation_error(stencil: list, interval: str = DEFAULT_INTERVAL):
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

    intp_eq = dot_product(coef, f_te)
    # calculate weighted sum of Taylor series.
    # for instance, 2nd-order 2-point interpolation is
    # intp_eq [= f(h)/2 + f(-h)/2)] = f(0) + f^(2)*h**2/2 + ...

    h = sp.symbols(interval)
    return sp.expand(
        sp.simplify(
            sp.symbols(DEFAULT_DIFFERENTIAND)
            - sp.nsimplify(intp_eq, rational=True, tolerance=1e-10)
        )
    ).as_leading_term(h)
    # extract the leading-order of errer term.
    # A interpolation formulation with error term is, for instance,
    # f(0) = (f(h) + f(-h))/2 - f^(2)*h**2/2 - ...
    # to extract error terms, reformulate intp_eq as
    # f(0) - intp_eq = - f^(2)*h**2/2 - ...
