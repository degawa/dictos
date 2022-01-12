def LagrangianBasis(x, degreeOfPolynomial, pointAt, xSet=None):
    import sympy as sp

    if degreeOfPolynomial <= 0:
        raise ValueError("degree of polynomial has to be greater than 0")

    numDataSet = degreeOfPolynomial+1

    if xSet is None:
        xSet = sp.symbols('x0:{:d}'.format(numDataSet))

    index = list(range(numDataSet))
    index.remove(pointAt)

    return sp.prod([(x-xSet[j])/(xSet[pointAt]-xSet[j]) for j in index])


def LagrangianPoly(x, xSet, fSet):
    if len(xSet) != len(fSet):
        raise ValueError(
            "The number of elements of xSet({:d}) and fSet({:d}) does not match."
            .format(len(xSet), len(fSet)))

    numDataSet = len(xSet)
    degreeOfPolynomial = numDataSet-1

    return sum([LagrangianBasis(x, degreeOfPolynomial, i, xSet)*fSet[i]
                for i in range(numDataSet)])


def Derivative(function, x, orderOfDifference=1):
    import sympy as sp

    return sp.simplify(
        sp.diff(function, x, orderOfDifference).subs([(x, 0)]))
