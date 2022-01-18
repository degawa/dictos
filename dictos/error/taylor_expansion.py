"""
Custom exceptions for errors related to Taylor Expansion.
"""


class TaylroExpansionError(Exception):
    """
    Base class for error related to Taylor Expansion.
    """

    pass


class UnsupportedOrderOfDerivativeError(TaylroExpansionError):
    """
    Exception raised for errors
    that specifying unsupported order of derivative.

    Attributes:
        deriv (int): order of derivative which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, deriv: int) -> None:
        self.message = (
            "The Taylor Expansion does not support "
            + "the derivative of the specified order {}. ".format(deriv)
            + "Specify an order of derivative greater than 0."
        )

    def __str__(self) -> str:
        return self.message