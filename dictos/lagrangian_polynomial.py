import sympy as sp


def lagrangian_basis(x, degree, point_at, x_set=None):
    """
    create Lagrangian basis from give coordinate symbols.

    Args:
        x (sympy symbol): symbol representing independent variable.
        degree (int): degree of polynomial.
        point_at (int): a index indicating the point
            where Lagrangian basis polynomial is defined in x_set.
        x_set (list or tuple of sympy symbols, optional):
            set of coordinate values. Defaults to None.

    Raises:
        ValueError: if degree of polynomial is less than 1.

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
    """calculate symbolically a lagrangian interpolation polynomial
    from ginve set of coordinates and functions.

    Args:
        x (sympy symbol): symbol representing independent variable.
        x_set (list or tuple of sympy symbols):
            set of coordinate values.
        f_set (list or tuple of sympy symbols):
            set of coordinate values.

    Raises:
        ValueError: if length of x_set and f_set are different.

    Returns:
        sympy Expr: a Lagrangian polynomial.

    Examples:
        >>> from dictos import lagrangian_polynomial as lp
        >>> import sympy as sp
        >>> x = sp.symbols("x")
        >>> dx = sp.symbols("dx")
        >>> x_set = [-dx, 0, dx]
        >>> f_set = sp.symbols("f0:{:d}".format(len(x_set)))
        >>> lp.lagrangian_poly(x, x_set, f_set)
        f0*x*(-dx + x)/(2*dx**2) - f1*(-dx + x)*(dx + x)/dx**2 + f2*x*(dx + x)/(2*dx**2)
    """
    if len(x_set) != len(f_set):
        raise ValueError(
            "The number of elements of xSet({:d}) and fSet({:d}) are different.".format(
                len(x_set), len(f_set)
            )
        )

    num_set = len(x_set)
    degree = num_set - 1
    # define an degree of polynomial.
    # (n-1)-th degree polynomial can be constrcuted using n points.

    return sum(
        [lagrangian_basis(x, degree, i, x_set) * f_set[i] for i in range(num_set)]
    )
    # calculate the linear combination of set of functions
    # at coordinates and the Lagrangian basis polynomials.


def derivative(expr, x, deriv=1):
    """calculate symbolically a derivative at x=0.

    Args:
        expr (sympy Expr): a derivative function.
        x (sympy symbol): symbol representing independent variable.
        deriv (int, optional): order of derivative. Defaults to 1.

    Returns:
        sympy Expr: a function differentiated symbolically at x=0.

    Examples:
        >>> from dictos import lagrangian_polynomial as lp
        >>> import sympy as sp
        >>> x, a, b, c = sp.symbols("x a b c")
        >>> lp.derivative(a*x**2+b*x+c, x, 1)
        b
        >>> lp.derivative(a*x**2+b*x+c, x, 2)
        2*a
    """
    return sp.simplify(sp.diff(expr, x, deriv).subs([(x, 0)]))
    # substitute 0 after differentiation
