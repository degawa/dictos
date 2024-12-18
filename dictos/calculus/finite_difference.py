import sympy as sp

from dictos.defaults import (
    DEFAULT_DIFFERENTIAND,
    DEFAULT_INTERVAL,
    DEFAULT_INDEPENDENT_VARIABLE,
)
from dictos.utilities.spec import (
    is_not_natural_number,
    is_valid_accuracy_order_for_generating_central_form,
)
from dictos.discrete.stencil import (
    create_coordinate_symbols,
    create_differentiand_symbols,
)
from dictos.utilities.utils import (
    simplify_coefficients,
    extract_coefficients_as_numer_denom,
    drop_coefficient_of_1,
    sort_by_subscript,
    is_odd,
    is_even,
)
from dictos.linalg.linalg import dot_product, div
from dictos.poly.lagrangian_polynomial import lagrangian_poly, derivative
from dictos.series.taylor_expansion import taylor_series, derivative_symbol
from dictos.calculus.exceptions import (
    UnsupportedOrderOfDerivativeError,
    InvalidOrderOfAccuracyForCentralFormError,
)
from dictos.core.expr import Expr
from dictos.core.grid_type import GridType


def equation(
    stencil: list,
    deriv: int = 1,
    interval: str = DEFAULT_INTERVAL,
    sort: bool = True,
    keep_zero: bool = False,
) -> Expr:
    """
    derive finite difference equation based on given stencil.

    Args:
        stencil (list of int): relative point numbers
            used for discretization.
        deriv (int, optional): order of derivative. Defaults to 1.
        interval (str, optional): an interval symbol like `dx`.
            Defaults to DEFAULT_INTERVAL.
        keep_zero (bool, optional):
            flag to keep terms multiplied by 0 in the result.
            Defaults to False.

    Returns:
        dictos Expr: derived finite difference equation.

    Examples:
        >>> from dictos import finite_difference as fd
        >>> fd.equation([-2, -1, 0, 1, 2], deriv=1)
        (f_{-2} - 8*f_{-1} + 8*f_{1} - f_{2})/(12*h)
        >>> fd.equation([-1.5, -0.5, 0, 0.5, 1.5], deriv=1)
        (f_{-1.5} - 27*f_{-0.5} + 27*f_{0.5} - f_{1.5})/(24*h)
    """

    x_set = create_coordinate_symbols(stencil, interval)
    f_set = create_differentiand_symbols(x_set, DEFAULT_DIFFERENTIAND)
    # create set of coordinate and differentiand symbols from stencil.
    # [-2, -1, 0, 1, 2] -> [-2*h, -h, 0, h, 2*h]
    # [-2*h, -h, 0, h, 2*h] -> [f_{-2}, f_{-1}, f_{0}, f_{1}, f_{2}]

    # derive finite difference coefficients based on given stencil,
    # and then calculate dot product of the coefs and differentiands.
    if keep_zero:
        numer, denom = coefficients(stencil, deriv, as_numer_denom=True)
        eq = div(
            drop_coefficient_of_1(dot_product(numer, f_set, evaluate=False)),
            denom * sp.symbols(interval) ** deriv,
        )
        # get coefficients as numerator and denominator and then
        # calculate dot product of numerator and differentiands
        # without evaluation to keep the term with a coefficient of 0
        # and avoid reducing the coefficients.
        # `div` calculates division with evaluation,
        # but the coefficients are not reduced.
    else:
        coef = coefficients(stencil, deriv)
        eq = sp.simplify(dot_product(coef, f_set) / sp.symbols(interval) ** deriv)
        # The result of dot product is divided by interval symbol,
        # because `coef` does not contain interval symbol like `dx`.
        # Terms multiplied by a coefficient 0 is eliminated from a result like
        # (-8*f_{-1} + f_{-2} + 8*f_{1} - f_{2})/(12*h).
        # `keep_zero=True` keep terms  multiplied by a coefficient 0,
        # (f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h).

    if sort:
        eq = sort_by_subscript(eq)
    # sort numerator by subscript.

    return Expr(eq)


def coefficients(stencil: list, deriv: int = 1, as_numer_denom: bool = False):
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
    if is_not_natural_number(deriv):
        raise UnsupportedOrderOfDerivativeError(deriv)
        # raise error
        # - if unsupported order of derivative (deriv < 1)

    x_set = create_coordinate_symbols(stencil, DEFAULT_INTERVAL)
    f_set = create_differentiand_symbols(x_set, DEFAULT_DIFFERENTIAND)
    # create set of coordinate and differentiand symbols from stencil.
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


def truncation_error(stencil: list, deriv: int, interval: str = DEFAULT_INTERVAL):
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

    fd_eq = dot_product(coef, f_ts)
    # calculate weighted sum of Taylor series.
    # for instance, 2nd-order 3-point central finite difference
    # for 1st derivative is
    # fd_eq [= f(h)/2 - f(-h)/2)] = f^(1)*h + f^(3)*h**3/6 + ...

    h = sp.symbols(interval)
    return sp.expand(
        sp.simplify(
            derivative_symbol(DEFAULT_DIFFERENTIAND, deriv)
            - sp.nsimplify(fd_eq / h**deriv, rational=True, tolerance=1e-10)
        )
    ).as_leading_term(h)
    # extract the leading-order of errer term.
    # A finite difference formulation with error term is, for instance,
    # f^(1) = (f(h) - f(-h))/(2*h) - f^(3)*h**3/6 - ...
    # to extract error terms, reformulate fd_eq as
    # f^(1) - fd_eq/h**1 = - f^(3)*h**3/6 - ...


def generate(
    deriv: int = 1,
    acc: int = 2,
    grid_type: GridType = GridType.REGULAR,
    as_equation: bool = False,
    consistent: bool = False,  # Reserved for future extension
):
    """
    generate a finite difference equation or a coefficient
    based on the order of derivative and the order of accuracy.

    Args:
        deriv (int, optional): Order of derivative. Defaults to 1.
        acc (int, optional): Order of accuracy. Must be even and >=2.
            Defaults to 2.
        grid_type (GridType, optional): Type of target grid system
            for generating equations or coefficients.
            Must be one of GridType.REGULAR, GridType.CELL_CENTERED, or GridType.STAGGERED.
            Defaults to GridType.REGULAR.
        as_equation (bool, optional): If True, returns equation in dictos Expr
            instead of coefficients in Tuple of sympy.Expr.
            Defaults to False.
        consistent (bool, optional): Reserved for futuer extension.
            Defaults to False.

    Returns:
        Union[dictos.Expr, Tuple[sympy.Expr]]: generated finite difference equation
        or generated coefficinets depending on `as_equation`.

    Raises:
        UnsupportedOrderOfDerivativeError: If deriv < 1
        ValueError: If acc is not an even number >= 2 or if grid_type is invalid

    Note:
        - If consistent is False, for staggered grid, uses regular grid method with even order derivatives
        and cell-centered grid method with odd order derivatives.
    """

    # validate order of derivative
    if is_not_natural_number(deriv):
        raise UnsupportedOrderOfDerivativeError(deriv)
        # raise error
        # - if unsupported order of derivative (deriv < 1)

    # validate order of accuracy
    if not is_valid_accuracy_order_for_generating_central_form(acc):
        raise InvalidOrderOfAccuracyForCentralFormError(acc)
        # raise error
        # - if acc is not positive and even

    if grid_type == GridType.REGULAR:
        return _generate_on_regular_grid(deriv, acc, as_equation)

    elif grid_type == GridType.CELL_CENTERED:
        return _generate_on_cell_centered_grid(deriv, acc, as_equation)

    elif grid_type == GridType.STAGGERED:
        return _generate_on_staggered_grid(deriv, acc, as_equation, consistent)

    else:
        raise ValueError(
            f"unsupported grid type: {grid_type}"
            "Must be one of: REGULAR, CELL_CENTERED, STAGGERED"
        )


def _generate_on_regular_grid(deriv: int, acc: int, as_equation: bool = False):
    """
    generate a finite difference equation or coefficients for regular grid system.

    Args:
        deriv (int): Order of derivative.
        acc (int): Order of accuracy.
        as_equation (bool, optional): If True, returns equation in dictos Expr
            instead of coefficients in Tuple of sympy.Expr.
            Defaults to False.

    Returns:
        Union[dictos.Expr, Tuple[sympy.Expr]]: generated finite difference equation
            or generated coefficinets depending on `as_equation`.

    Note:
        - Half width of the stencil is (deriv + acc - 1)//2 accoding to the table below:
            | acc | deriv | stencil width |
            |:---:|:-----:|:-------------:|
            |  2  |   1   |       3       |
            |  2  |   2   |       3       |
            |  2  |   3   |       5       |
            |  2  |   4   |       5       |
            |  2  |   5   |       7       |
            |  4  |   1   |       5       |
            |  4  |   2   |       5       |
            |  4  |   3   |       7       |
            |  4  |   4   |       7       |
            |  4  |   5   |       9       |
        - The stencil width is 2 * half_width + 1, where +1 accounts for the center point.
    """

    stencil_half_width = (deriv + acc - 1) // 2
    # calculate half-width of the stencil

    stencil = [i for i in range(-stencil_half_width, stencil_half_width + 1)]
    # generate stencil [-half_width, ..., 0, ..., half_width]
    # stencil width is `stencil_half_width * 2 + 1`
    # where +1 is for point 0

    if as_equation:
        return equation(stencil, deriv)
    else:
        return coefficients(stencil, deriv)


def _generate_on_cell_centered_grid(deriv: int, acc: int, as_equation: bool = False):
    """
    generate a finite difference equation or coefficients for cell-centered grid system.

    Args:
        deriv (int): Order of derivative.
        acc (int): Order of accuracy.
        as_equation (bool, optional): If True, returns equation in dictos Expr
            instead of coefficients in Tuple of sympy.Expr.
            Defaults to False.

    Returns:
        Union[dictos.Expr, Tuple[sympy.Expr]]: generated finite difference equation
            or generated coefficinets depending on `as_equation`.

    Note:
        - Half width of the stencil  is (deriv + acc)//2 accoding to the table below:
            | acc | deriv | stencil width |
            |:---:|:-----:|:-------------:|
            |  2  |   1   |       2       |
            |  2  |   2   |       4       |
            |  2  |   3   |       4       |
            |  2  |   4   |       6       |
            |  2  |   5   |       6       |
            |  4  |   1   |       4       |
            |  4  |   2   |       6       |
            |  4  |   3   |       6       |
            |  4  |   4   |       8       |
            |  4  |   5   |       8       |
        - The stencil width is 2 * half_width.
    """

    stencil_half_width = (deriv + acc) // 2
    # calculate half-width of the stencil

    stencil = [i + 0.5 for i in range(-stencil_half_width, stencil_half_width)]
    # generate stencil [-half_width+0.5, ..., -0.5, 0.5, ..., half_width-0.5]
    # stencil width is `stencil_half_width * 2`

    if as_equation:
        return equation(stencil, deriv)
    else:
        return coefficients(stencil, deriv)


def _generate_on_staggered_grid(
    deriv, acc, as_equation: bool = False, consistent: bool = False
):
    """
    generate a finite difference equation or coefficients for staggered grid system.

    Args:
        deriv (int): Order of derivative.
        acc (int): Order of accuracy.
        as_equation (bool, optional): If True, returns equation in dictos Expr
            instead of coefficients in Tuple of sympy.Expr.
            Defaults to False.
        consistent (bool, optional): Reserved for future extension. Defaults to False.

    Returns:
        Union[dictos.Expr, Tuple[sympy.Expr]]: generated finite difference equation
            or generated coefficinets depending on `as_equation`.

    Note:
        - For even order derivatives, uses regular grid method,
        and for odd order derivatives, uses cell-centered grid method.
        See the table below:
            | acc | deriv | stencil width | def. point |
            |:---:|:-----:|:-------------:|:----------:|
            |  2  |   1   |       2       |     mid    |
            |  2  |   2   |       3       |    grid    |
            |  2  |   3   |       4       |     mid    |
            |  2  |   4   |       5       |    grid    |
            |  2  |   5   |       6       |     mid    |
            |  4  |   1   |       4       |     mid    |
            |  4  |   2   |       5       |    grid    |
            |  4  |   3   |       6       |     mid    |
            |  4  |   4   |       7       |    grid    |
            |  4  |   5   |       8       |     mid    |
    """

    if is_even(deriv):
        return _generate_on_regular_grid(deriv, acc, as_equation)
    else:
        return _generate_on_cell_centered_grid(deriv, acc, as_equation)
