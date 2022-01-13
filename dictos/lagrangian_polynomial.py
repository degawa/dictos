import sympy as sp


def lagrangian_basis(x, degree, point_at, x_set=None):
    """
    create Lagrangian basis from give coordinate symbols.

    Args:
        x (sympy symbol): symbol representing independent variable
        degree (int): degree of polynomial
        point_at (int): a index indicating the point
            where Lagrangian basis polynomial is defined in x_set.
        x_set (list or tuple of sympy symbols, optional):
            set of coordinate values. Defaults to None.

    Raises:
        ValueError: if degree of polynomial is less than 1

    Returns:
        sympy Expr: a Lagrange basis polynomial $l(x)|_{xset[point_at]}$.

    Examples:
        >>> from dictos import lagrangian_polynomial as lp
        >>> import sympy as sp
        >>> x = sp.symbols("x")
        >>> lp.lagrangian_basis(x, degree=3, point_at=0)
        (x - x1)*(x - x2)*(x - x3)/((x0 - x1)*(x0 - x2)*(x0 - x3))
    """
    if degree <= 0:
        raise ValueError("degree of polynomial has to be greater than 0")

    # TODO: raise error when at least one coordinate values in the x_set appears more than once.

    num_set = degree + 1
    # n+1 points are required to construct an n-th degree polynomial.

    if x_set is None:
        x_set = sp.symbols("x0:{:d}".format(num_set))
        # create set of coordinate values like `(x0, x1, x2, ...)`

    index = list(range(num_set))
    index.remove(point_at)
    # create indices indicating set of coordinate values
    # used to construct the polynomial
    # once create `[0, 1, 2, ...]`
    # and then remove an index, `0`` for example, indicating the point
    # where Lagrangian basis polynomial is defined

    return sp.prod([(x - x_set[j]) / (x_set[point_at] - x_set[j]) for j in index])
    # calculate a Lagrangian basis polynomial.
    # `j` is used as an index. When constructing a Lagrangian polynomial,
    # a Lagrangian basis polynomial is constructed at several points.
    # `i` is used as an index to indicate those points, so `j` is used here.


def lagrangian_poly(x, x_set, f_set):
    if len(x_set) != len(f_set):
        raise ValueError(
            "The number of elements of xSet({:d}) and fSet({:d}) are different.".format(
                len(x_set), len(f_set)
            )
        )

    num_set = len(x_set)
    degree = num_set - 1

    return sum(
        [lagrangian_basis(x, degree, i, x_set) * f_set[i] for i in range(num_set)]
    )


def Derivative(function, x, orderOfDifference=1):
    return sp.simplify(sp.diff(function, x, orderOfDifference).subs([(x, 0)]))
