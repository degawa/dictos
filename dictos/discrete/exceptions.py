"""
Custom exceptions for errors related to stencil.
"""

from dictos.utilities.exceptions.utils import find_duplicated_points


class StencilError(Exception):
    """
    Base class for error related to stencil.
    """

    pass


class TooNarrowError(StencilError):
    """
    Exception raised for errors that the stencil is too narrow.
    """

    def __init__(self, stencil) -> None:
        self.message = (
            f"The stencil width {len(stencil)} is too narrow. "
            + "Use a stencil containing at least 2 points."
        )

    def __str__(self) -> str:
        return self.message


class DuplicatedPointError(StencilError):
    """
    Exception raised for errors
    that at least a point is duplicated in the stencil.

    Attributes:
        stencil (list of int or list of float):
            stencil which caused the error.
        message (str): Explanation of the error.
    """

    def __init__(self, stencil) -> None:
        self.stencil = stencil
        self.message = (
            "The stencil has at least one duplicated point. "
            + "Remove/fix duplicated point(s) "
            + f"{find_duplicated_points(self.stencil)} from the stencil."
        )

    def __str__(self) -> str:
        return self.message


class ContainsZeroError(StencilError):
    """
    Exception raised for errors
    that the stencil contains `0`.

    Attributes:
        message (str): Explanation of the error.
    """

    def __init__(self) -> None:
        self.message = (
            "The stencil contains `0`. "
            "Remove `0` from the stencil using `.remove(0)` etc."
        )

    def __str__(self) -> str:
        return self.message
