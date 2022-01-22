import sympy as sp

from .spec import (
    are_different_length,
)
from .error.linear_algebra import InconsistentDataSetError


def dot_product(vec1, vec2, evaluate=True):
    """
    calculate dot product of two list
    containing sympy numbers and symbols without evaluation.

    Args:
        vec1 (list of int or float, or sympy Mul): 1D-list.
            The list length must be equal to it of vec2.
        vec2 (list of int or float, or sympy Mul): 1D-list.
            The list length must be equal to it of vec1.
        evaluate (bool, optional): flag to evaluate the result.
            Defaults to True.

    Raises:
        InconsistentDataSetError: if two lists are inconsistent.

    Returns:
        sympy Expr: dot product of the passed two lists.
    """
    if are_different_length(vec1, vec2):
        raise InconsistentDataSetError(vec1, vec2)
        # raise error if
        # two lists are inconsistent.

    begin_ = len(vec2) - 1  # exclude first term
    end_ = -1  # to generate numbers up to 0 using range()
    step_ = -1
    eq = sp.Mul(vec1[begin_], vec2[begin_], evaluate=evaluate)
    for i in range(begin_ + step_, end_, step_):
        eq = sp.Add(eq, sp.Mul(vec1[i], vec2[i], evaluate=evaluate), evaluate=evaluate)
    # there is no way to accumulate a list of sympy symbols
    # without evalulation.
    # So for-loop and sympy.Add and .Mul are used.
    # The for-loop is downstepped
    # so that the result are sorted when it is printed.
    # The result of the accumulation, `((-a*f_{-1} +b*f_{0}) + a*f_{1})`,
    # is printed as `a*f_{1} + b*f_{0} - a*f_{-1}`.
    # downstepping is used for sorting.

    return eq


def div(numer, denom):
    """
    calculate division of two sympy Expr.
    The result is always evaluated.

    Args:
        numer (sympy Expr): Algebraic expression to be devided
        denom (sympy Expr): Algebraic expression to divide.
            The Expr must be a single term.

    Returns:
        sympy Expr: result of division eq/denom.
    """
    return sp.Mul(numer, 1 / denom)
    # calculate Mul with evaluate=True,
    # because the result with evaluate=False will be like
    # `(1/(12*h))*(f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})`,
    # unlike the result we want as follows:
    # `(f_{-2} - 8*f_{-1} + 0*f_{0} + 8*f_{1} - f_{2})/(12*h)`
