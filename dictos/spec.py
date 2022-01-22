"""
Describe specification of dictos and provide functions for inspections.
"""

_MINIMUM_STENCIL_WIDTH = 2


def narrower_than_minimum_width(data_set) -> bool:
    """
    Returns True
    if stencil used in data_set is narrower than
    the minimum value specified in the spec.py.

    Args:
        data_set (list): data set to be checked.

    Returns:
        bool: True if stencil used in data_set is narrower.
    """
    return len(data_set) < _MINIMUM_STENCIL_WIDTH


def is_not_natrual_number(number, include_zero: bool = False) -> bool:
    """
    Returns True if `number` is not the natural number.

    Args:
        number : an integer number to be checked.
        include_zero (bool, optional): A flag to swith to include 0
            in the natural numbers.
            Defaults to False.

    Returns:
        bool: True if `number` is not the natural number.
    """
    return not is_natural_number(number, include_zero)


def has_zero(stencil) -> bool:
    """
    Return True if the stencil has 0.

    Args:
        stencil (list of int or float): stencil on regular or
            staggered grid.

    Returns:
        bool: True if the stencil has 0.
    """
    return 0 in stencil


def has_duplicated_points(stencil) -> bool:
    """
    Returns True if there is at least one duplicated points in the stencil.

    Args:
        stencil (list of int or float): stencil on regular or
            staggered grid.

    Returns:
        bool: True if there is at least one duplicated points in the stencil.
    """
    return len(stencil) != len(set(stencil))


def are_different_length(list1, list2) -> bool:
    """Return True if length of two lists are different.

    Args:
        list1 (list): a list to be checked.
        list2 (list): a list to be checked.

    Returns:
        bool: True if length of two lists are different.
    """
    return not are_same_length(list1, list2)


def is_not_assumed_length(list, assumed_length: int) -> bool:
    """
    Return Ture if the list length is not the assumed length.

    Args:
        list (list): a list to be checked.
        assumed_length (int): assumed length of the list.

    Returns:
        bool: True if the list length is not the assumed length.
    """
    return not is_assumed_length(list, assumed_length)


def is_natural_number(number, include_zero: bool = False) -> bool:
    """
    Returns True if `number` is the natural number.

    Args:
        number : an integer number to be checked.
        include_zero (bool, optional): A flag to swith to include 0
            in the natural numbers.
            Defaults to False.

    Returns:
        bool: True if `number` is the natural number.
    """
    lower_bound = 0 if include_zero else 1
    return isinstance(number, int) and number >= lower_bound


def are_same_length(list1, list2) -> bool:
    """Return True if length of two lists are the same.

    Args:
        list1 (list): a list to be checked.
        list2 (list): a list to be checked.

    Returns:
        bool: True if length of two lists are the same.
    """
    return len(list1) == len(list2)


def is_assumed_length(list, assumed_length: int) -> bool:
    """
    Return Ture if the list length is the assumed length.

    Args:
        list (list): a list to be checked.
        assumed_length (int): assumed length of the list.

    Returns:
        bool: True if the list length is the assumed length.
    """
    return len(list) == assumed_length
