import sympy as sp

from dictos.utilities.spec import are_different_length
from dictos.linalg.exceptions import InconsistentDataSetError


def dot_product(vec1, vec2, evaluate: bool = True):
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

    if evaluate:
        eq = sum([vec1[i] * vec2[i] for i in range(len(vec1))])
    else:
        terms = [sp.Mul(vec1[i], vec2[i], evaluate=False) for i in range(len(vec1))]
        eq = sp.Add(*terms, evaluate=False)
        # there is no way to accumulate a list of sympy symbols
        # without evaluation.
        # So for-loop and sympy.Add and .Mul are used.

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
