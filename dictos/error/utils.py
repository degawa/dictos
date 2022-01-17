import collections


def find_duplicated_points(stencil: list) -> list:
    """Find duplicated points from the stencil.

    Args:
        stencil (list of int or float): stencil which caused the error.

    Returns:
        list of int or float: duplicated points number.
    """
    return [point for point, count in collections.Counter(stencil).items() if count > 1]
