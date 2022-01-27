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
            + f"the derivative of the specified order {deriv}. "
            + "Specify an order of derivative greater than 0."
        )

    def __str__(self) -> str:
        return self.message


class NumberOfExpansionTermsIsNotNaturalNumberError(TaylroExpansionError):
    """
    Exception raised for errors
    that the number of series expansion terms is not the natural number.

    Attributes:
        term (int): number of terms which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, term: int) -> None:
        self.message = (
            f"The number of series expansion terms {term} "
            + "is not the natural number. "
            + "Specify the number of terms greater than 0."
        )

    def __str__(self) -> str:
        return self.message
