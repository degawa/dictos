import sympy as sp
from typing import List

from dictos.utilities.utils import simplify_coefficients, is_odd, sort_by_subscript
from dictos.calculus import finite_difference as fd
from dictos.core.grid_type import GridType
from dictos.discrete.stencil import (
    create_coordinate_symbols,
    create_differentiand_symbols,
)
from dictos.linalg.linalg import dot_product, add, scale
from dictos.core.expr import Expr


def generate(acc: int, as_numer_denom: bool = False, as_equation: bool = False):
    """
    generate the equation or coefficients
    for the linear filter on the regular grid

    Args:
        acc (int): Order of accuracy (must be even and positive)
        as_numer_denom (bool): If True, return coefficients as numerator/denominator
        as_equation (bool): If True, return as symbolic equation

    Returns:
        Union[sp.Expr, Expr]: equation or coefficient for linear filter
    """

    # validate order of accuracy
    if is_odd(acc) or acc <= 0:
        raise ValueError(
            f"order of accuracy `acc` must be an even number >=2, got: {acc}"
        )

    # Handle conflict flags
    if as_numer_denom and as_equation:
        as_equation = False

    # generate finite difference coefficients and calculate damping coefficients
    denom = 2**acc
    exponenet = (acc - 2) // 2

    damp = fd.generate(deriv=acc, acc=2, grid_type=GridType.REGULAR)
    damp = scale(scale(damp, sp.Pow(denom, -1)), (-1) ** exponenet)
    # generate 2nd-order *acc-th* derivative finite difference in central form
    # and then calculate (-1)**((acc-2)//2)*(h/2)**acc * ∂**acc f/∂h**acc |_i=0
    # h is cancelled by with finite difference form of ∂**acc f/∂h**acc

    stencil_width = len(damp)
    main_component = [0] * stencil_width  # create zero list [0, 0, 0, 0, 0]
    main_component[stencil_width // 2] = 1  # update list to [0, 0, 1, 0, 0] for f

    # combine components to construct filter
    coef = add(damp, main_component)
    # f_filtered = f + (-1)**((acc-2)//2)*(h/2)**acc * ∂**acc f/∂h**acc

    if as_equation:
        return _generate_equation(coef, stencil_width)
    else:
        return simplify_coefficients(coef, as_numer_denom=as_numer_denom)


def _generate_equation(coefficients: List[sp.Expr], stencil_width: int) -> Expr:
    """
    generate symbolic equation from coefficients

    Args:
        coefficients (List[sp.Expr]): List of coefficients
        stencil_width (int): Width of the stencil

    Returns:
        Expr: Symbolic equation
    """
    half_width = stencil_width // 2
    stencil_range = range(-half_width, half_width + 1)
    x_set = create_coordinate_symbols(list(stencil_range))
    f_set = create_differentiand_symbols(x_set)
    eq = sp.simplify(dot_product(coefficients, f_set))
    eq = sort_by_subscript(eq)

    return Expr(eq)
