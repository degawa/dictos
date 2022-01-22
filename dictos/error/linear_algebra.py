"""
Custom exceptions for errors related to Linear Algebra.
"""


class LinearAlgebraError(Exception):
    """
    Base class for error related to Linear Algebra.
    """

    pass


class InconsistentDataSetError(LinearAlgebraError):
    """
    Exception raised for errors
    that data set of two vectors are inconsistent.

    Attributes:
        vec1 (list of int or float, or sympy Mul): a vector
            which caused the error.
        vec2 (list of int or float, or sympy Mul): a vector
            which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, vec1: list, vec2: list) -> None:
        self.vec1 = vec1
        self.vec2 = vec2
        self.message = (
            "Inconsistent data set of two vectors. "
            + "The number of elements of "
            + "vector 1 [{:d}] and vector 2 [{:d}] are different.".format(
                len(self.vec1), len(self.vec2)
            )
        )

    def __str__(self) -> str:
        return self.message
