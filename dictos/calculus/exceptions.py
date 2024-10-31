"""
Custom exceptions for errors related to Finite Difference.
"""


class FiniteDifferenceError(Exception):
    """
    Base class for error related to Finite Difference.
    """

    pass


class UnsupportedOrderOfDerivativeError(FiniteDifferenceError):
    """
    Exception raised for errors
    that specifying unsupported order of derivative.

    Attributes:
        deriv (int): order of derivative which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, deriv: int) -> None:
        self.message = (
            "The Finite Difference Method does not support "
            + f"the derivative of the specified order {deriv}. "
            + "Specify an order of derivative greater than 1."
        )

    def __str__(self) -> str:
        return self.message
