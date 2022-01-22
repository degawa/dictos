"""
Custom exceptions for errors related to Lagrangian Polynomial.
"""


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

    def __init__(self, x_set, f_set) -> None:
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


class DegreeOfPolynomialIsNotNaturalNumberError(LagrangianPolynomialError):
    """
    Exception raised for errors
    that the degree of polynomial is not the naturarl number.

    Attributes:
        degree (int): degree of polynomial which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, degree) -> None:
        self.degree = degree
        self.message = (
            "The degree of polynomial ({}) is not the natural number. ".format(degree)
            + "Specify the degree of polynomial greater than 0."
        )

    def __str__(self) -> str:
        return self.message


class InconsistentDataSetAndDegreeOfPolynomialError(LagrangianPolynomialError):
    """
    Exception raised for errors
    that there is an inconsistency between the degree of polynomial
    and number of data set.

    Attributes:
        degree (int): degree of polynomial which causedthe error.
        x_set (list of int or float): set of independend variable
            which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, degree: int, x_set: list) -> None:
        self.degree = degree
        self.x_set = x_set
        self.message = (
            "Inconsistency between the degree of polynomial and data set. "
            + "The degree of polynomial ({}) must be the number of data set ({}) - 1.".format(
                degree, len(x_set)
            )
        )

    def __str__(self) -> str:
        return self.message
