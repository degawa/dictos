import sympy as sp
from typing import List, Callable
from operator import add as add_op, mul as mul_op
import sympy as sp

from dictos.core.expr import Expr
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


def scale(x, factor) -> List:
    """
    multiply each element in the list by a scalar value

    Args:
        x (List): input list
        factor (int or float): scalar multiplier

    Returns:
        List: Scaled list

    Examples:
        >>> scale([1, 2, 3], 2)
        [2, 4, 6]
    """
    return _scalar_op(mul_op, x, factor)


def add(a: List, b: List) -> List:
    """
    calculate element-wise sum of two lists

    Args:
        a (List): First input list
        b (List): Second input list (must be the same length as a)

    Returns:
        List: Element-wise sum

    Raises:
        ValueError: If input lists have different lengths

    Examples:
        >>> add([1, 2, 3], [4, 5, 6])
        [5, 7, 9]
    """
    if len(a) != len(b):
        raise ValueError(f"Lists must have same length, got {len(a)} and {len(b)}")

    return _list_op(add_op, a, b)


def _list_op(op: Callable, *lists: List) -> List:
    """
    apply element-wise operation to multiple lists

    Args:
        op (Callable): Binary operation function
        *lists: Variable number of input lists

    Returns:
        List: List containing operation results

    Raises:
        ValueError: If input lists have different lengths

    Examples:
        >>> list_op(add_op, [1, 2], [3, 4], [5, 6])
        [9, 12]
    """
    if len(set(len(lst) for lst in lists)) > 1:
        raise ValueError("All input lists must have the same length")

    return list(map(op, *lists))


def _scalar_op(op: Callable, lst: List, scalar) -> List:
    """
    apply operation between each list element and a scalar value

    Args:
        op (Callable): Binary operation function
        lst (List): Input list
        scalar (int or float): Scalar value

    Returns:
        List: List containing operation results

    Examples:
        >>> scalar_op(mul_op, [1, 2, 3], 2)
        [2, 4, 6]
    """
    return list(map(lambda x: op(x, scalar), lst))
