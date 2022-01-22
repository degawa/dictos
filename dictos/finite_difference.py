import sympy as sp

from .spec import is_not_natrual_number
from .utils import (
    DEFAULT_FUNCTION,
    DEFAULT_INTERVAL,
    DEFAULT_INDEPENDENT_VARIABLE,
    create_coordinate_symbols,
    create_function_symbols,
    simplify_coefficients,
    extract_coefficients_as_numer_denom,
)
from .linalg import dot_product, div
from .lagrangian_polynomial import lagrangian_poly, derivative
from .taylor_expansion import taylor_series, derivative_symbol
from .error.finite_difference import UnsupportedOrderOfDerivativeError


def equation(
    stencil,
    deriv=1,
    interval=DEFAULT_INTERVAL,
    same_subscripts_as_stencil=False,
    evaluate=True,
):
    """
    derive finite difference equation based on given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        deriv (int, optional): order of derivative. Defaults to 1.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL.
        same_subscripts_as_stencil (bool, optional): flag
            to make function subscripts the same as the stencil.
            Defaults to False, the subscripts start from 0.
        evaluate (bool, optional): flag to evaluate the result.
            Defaults to True.

    Returns:
        sympy Expr: derived finite difference equation.

    Examples:
        >>> from dictos import finite_difference as fd
        >>> fd.equation([-2, -1, 0, 1, 2], deriv=1)
        (f_0 - 8*f_1 + 8*f_3 - f_4)/(12*h)
        >>> fd.equation([-1.5, -0.5, 0, 0.5, 1.5], deriv=1)
        (f_0 - 27*f_1 + 27*f_3 - f_4)/(24*h)
    """

    x_set = create_coordinate_symbols(stencil, interval)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION, same_subscripts_as_stencil)
    # create set of coordinate and function symbols from stencil.
    # [-2, -1, 0, 1, 2] -> [-2*h, -h, 0, h, 2*h]
    # [-2*h, -h, 0, h, 2*h] -> [f_0, f_1, f_2, f_3, f_4]
    #                          [f_{-2}, f_{-1}, f_{0}, f_{1}, f_{2}]
    #                          (same_subscripts_as_stencil = True)

    # derive finite difference coefficients based on given stencil,
    # and then calculate dot product of the coef and function values.
    if evaluate:
        coef = coefficients(stencil, deriv)
        eq = sp.simplify(
            sum([coef[i] * f_set[i] for i in range(len(f_set))])
            / sp.symbols(interval) ** deriv
        )
        # The result of dot product is divided by interval symbol,
        # because coef does not contain interval symbol like `dx`.
        # In many case, the results obtained by above operations is
        # what you may want like (f0 - 8*f1 + 8*f3 - f4)/(12*dx),
        # but in rare cases, an equation cannnot be expressed
        # by a rational expression,
        # (0.0416666666666667*f0 - 1.125*f1 + 1.125*f3 - 0.0416666666666667*f4)/h.
        # In such cases, specify `False` to the argument `evaluate`.
        # Another case, when same_subscripts_as_stencil is set `True`,
        # a resultant equation is not sorted like
        # (-8*f_{-1} + f_{-2} + 8*f_{1} - f_{2})/(12*h).
        # `evaluate=False` derives an sorted expression you may want,
        # (f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h).
    else:
        numer, denom = coefficients(stencil, deriv, as_numer_denom=True)
        eq = div(
            dot_product(numer, f_set, evaluate=False),
            denom * sp.symbols(interval) ** deriv,
        )
        # get coefficients as numerator and denominator and then
        # calculate dot product of numerator and function values
        # without evaluation to keep the term with a coefficient of 0
        # and avoid reducing the coefficients.
        # `div` calculates division with evaluation,
        # but the coefficients are not reduced.

    return eq


def coefficients(stencil, deriv=1, as_numer_denom=False):
    """
    derive finite difference coefficients based on given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        deriv (int, optional): order of derivative. Defaults to 1.
        as_numer_denom (bool, optional): flag to return the numerator
            and denominator separately.
            Defaults to False.

    Raises:
        UnsupportedOrderOfDerivativeError: if
            unsupported order of derivative (deriv < 1) is passed.

    Returns:
        list of sympy Rational: simplified coefficients.
            or
        list of sympy numbers, int:
            numerator and denominator of coefficients.
            coefficients are commutative
            with the least common multiple of the denominator.

    Examples:
        >>> from dictos import finite_difference as fd
        >>> fd.coefficients([-2, -1, 0, 1, 2], deriv=1)
        [1/12, -2/3, 0, 2/3, -1/12]
        >>> fd.coefficients([-1.5, -0.5, 0, 0.5, 1.5], deriv=1)
        [1/24, -9/8, 0, 9/8, -1/24]
        >>> fd.coefficients([-1.5, -0.5, 0, 0.5, 1.5], deriv=1, as_numer_denom=True)
        ([1, -27, 0, 27, -1], 24)
    """
    if is_not_natrual_number(deriv):
        raise UnsupportedOrderOfDerivativeError(deriv)
        # raise error
        # - if unsupported order of derivative (deriv < 1)

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_function_symbols(x_set, DEFAULT_FUNCTION)
    # create set of coordinate and function symbols from stencil.
    # [-2, -1, 0, 1, 2] -> [-2*h, -h, 0, h, 2*h]
    # [-2*h, -h, 0, h, 2*h] -> [f_0, f_1, f_2, f_3, f_4]

    x = sp.symbols(DEFAULT_INDEPENDENT_VARIABLE)
    poly = lagrangian_poly(x, x_set, f_set)
    numer_coef, denom_coef = extract_coefficients_as_numer_denom(poly, f_set)
    # extract numerator and denomitaor from the polynomial

    coef = [derivative(num / denom_coef[0], x, deriv) for num in numer_coef]
    # get coefficients of each terms as a list. Another expression
    # `derivative(lagrangian_poly(x, x_set, f_set), x, deriv).as_poly(f_set).coeffs()`
    # erases terms with a coefficient of 0.

    return simplify_coefficients(coef, as_numer_denom)
    # simplify floating-point number coefficients to ratioanl numbers


def truncation_error(stencil, deriv, interval=DEFAULT_INTERVAL):
    """
    derive the leading-order of error term
    in the finite difference equation based on the given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        deriv (int): order of derivative.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL.

    Returns:
        sympy Expr: the leading-order of error term

    Examples:
        >>> from dictos import finite_difference as fd
        >>> fd.truncation_error([-1, 0, 1], deriv=1)
        -f^(3)*h**2/6
        >>> fd.truncation_error([-1, 0, 1], deriv=2)
        -f^(4)*h**2/12
        >>> fd.truncation_error([-2, -1, 0, 1, 2], deriv=1)
        f^(5)*h**4/30
        >>> fd.truncation_error([-2, -1, 0, 1, 2], deriv=2)
        f^(6)*h**4/90
    """

    coef = coefficients(stencil, deriv)
    # derive finite difference coefficients based on given stencil

    x_set = create_coordinate_symbols(stencil, interval=interval)
    # create set of coordinate symbols from stencil.
    # [-2, -1, 0, 1, 2] -> [-2*h, -h, 0, h, 2*h]

    num_term = len(x_set) + deriv
    f_ts = [taylor_series(x, num_term) for x in x_set]
    # calculate Taylor series around points in x_set.

    fd_eq = sum([coef[i] * f_ts[i] for i in range(len(x_set))])
    # calculate weighted sum of Taylor series.
    # for instance, 2nd-order 3-point central finite difference
    # for 1st derivative is
    # fd_eq [= f(h)/2 - f(-h)/2)] = f^(1)*h + f^(3)*h**3/6 + ...

    h = sp.symbols(interval)
    return sp.expand(
        sp.simplify(
            derivative_symbol(DEFAULT_FUNCTION, deriv)
            - sp.nsimplify(fd_eq / h ** deriv, rational=True, tolerance=1e-10)
        )
    ).as_leading_term(h)
    # extract the leading-order of errer term.
    # A finite difference formulation with error term is, for instance,
    # f^(1) = (f(h) - f(-h))/(2*h) - f^(3)*h**3/6 - ...
    # to extract error terms, reformulate fd_eq as
    # f^(1) - fd_eq/h**1 = - f^(3)*h**3/6 - ...
