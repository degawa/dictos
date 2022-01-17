import collections


def extract_dupulicated_points(stencil: list) -> list:
    """Extract duplicated points from the stencil.

    Args:
        stencil (list of int or float): stencil which caused the error.

    Returns:
        list of int or float: dupulicated points number.
    """
    return [point for point, count in collections.Counter(stencil).items() if count > 1]
