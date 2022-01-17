class LagrangianPolynomialError(Exception):
    """
    Base class for error related to Lagrangian Polynomial.
    """

    pass


class InconsistentDataSetError(LagrangianPolynomialError):
    """
    Exception raised for errors
    that data set of independent variable and function are inconsistent.

    Attributes:
        x_set (list of int or float): set of independend variable
            which caused the error.
        f_set (list of int or float): set of function which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, x_set: list, f_set: list) -> None:
        self.x_set = x_set
        self.f_set = f_set
        self.message = (
            "Inconsistent data set of independent variable and function "
            + "The number of elements of x_set[{:d}] and f_fet[{:d}] are different.".format(
                len(x_set), len(f_set)
            )
        )

    def __str__(self) -> str:
        return self.message


class ViolateDegreeOfPolynomialAssumption(LagrangianPolynomialError):
    """
    Exception raised for errors
    that violate a theoretical assumption imposed
    on the degree of polynomial (degree > 0)

    Attributes:
        degree (int): degree of polynomial which causedthe error.
        message (str): Explanation of the error.
    """

    def __init__(self, degree: int) -> None:
        self.degree = degree
        self.message = (
            "Violates the theoretical assumption "
            + "imposed on the degree of polynomial. "
            + "The degree of polynomial has to be greater than 0 (passed value: {})".format(
                degree
            )
        )

    def __str__(self) -> str:
        return self.message
