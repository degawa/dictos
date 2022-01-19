"""
Custom exceptions for errors which caused unexpected results returned from Python.
The term "Internal" is defined for features implemented in Python and libraries.
"""


class InternalError(Exception):
    """
    Base class for error related to results from Python.
    """

    pass


class UnexpectedDenominatorError(InternalError):
    """
    Exception raised for errors
    that a denominator of an unexpected data type is returned from sympy.

    Attributes:
        denom (list of int of float): denominator which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, denom) -> None:
        length = ""
        if type(denom) is list or type(denom) is tuple:
            length = "[{}]".format(len(denom))
        self.message = (
            "Unexpected data type of a denominator "
            + str(type(denom))
            + length
            + ". "
            + "Confirm input/output data to/from sympy function(s) "
            + "and its internal operation."
        )

    def __str__(self) -> str:
        return self.message
