import sympy as sp


def lagrangian_basis(x, degree, point_at, x_set=None):
    if degree <= 0:
        raise ValueError("degree of polynomial has to be greater than 0")

    num_set = degree + 1

    if x_set is None:
        x_set = sp.symbols("x0:{:d}".format(num_set))

    index = list(range(num_set))
    index.remove(point_at)

    return sp.prod([(x - x_set[j]) / (x_set[point_at] - x_set[j]) for j in index])


def LagrangianPoly(x, xSet, fSet):
    if len(xSet) != len(fSet):
        raise ValueError(
            "The number of elements of xSet({:d}) and fSet({:d}) does not match.".format(
                len(xSet), len(fSet)
            )
        )

    numDataSet = len(xSet)
    degreeOfPolynomial = numDataSet - 1

    return sum(
        [
            lagrangian_basis(x, degreeOfPolynomial, i, xSet) * fSet[i]
            for i in range(numDataSet)
        ]
    )


def Derivative(function, x, orderOfDifference=1):
    return sp.simplify(sp.diff(function, x, orderOfDifference).subs([(x, 0)]))
