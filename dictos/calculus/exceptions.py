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


class InvalidOrderOfAccuracyForCentralFormError(FiniteDifferenceError):
    """
    Exception raised for errors
    that order of accucary is not even-order.

    Attributes:
        acc (int): order of accuracy which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, acc: int) -> None:
        self.message = (
            f"The order of accucary {acc} must be an even and positive number >=2. "
        )

    def __str__(self) -> str:
        return self.message
