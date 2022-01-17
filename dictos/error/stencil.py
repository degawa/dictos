class StencilError(Exception):
    """
    Base class for error related to stencil.
    """

    pass


class ContainsZeroError(StencilError):
    """
    Exception raised for errors
    that the stencil contains `0`.

    Attributes:
        message (str): Explanation of the error.
    """

    def __init__(self):
        self.message = (
            "The stencil contains `0`. "
            "Remove `0` from the stencil using `.remove(0)` etc."
        )

    def __str__(self):
        return self.message
