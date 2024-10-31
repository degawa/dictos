"""random input generator for tests
"""

import random
import string
from typing import Optional, List


STENCIL_HALF_WIDTH = 20  # up to 20th order accuracy
MAX_SYMBOL_LENGTH = 5


def random_string(len):
    """generate n-length random string

    Args:
        len (int): length of string
    """

    return "".join(random.choices(string.ascii_letters + string.digits, k=len))


def random_int(min: int, max: int, exclude: Optional[List[int]] = None) -> List[int]:
    """generate random list of integers from `min` to `max`.
    Integers listed in `exclude` are excluded from the generated list.

    Args:
        min (int): minimum value of integer in the list.
        max (int): maximum value of integer in the list.
        exclude (list of int): integers to be excluded from the list.
            Defaults to None.
    """
    num = list(range(min, max + 1))
    # +1 is the correction for exclusive stop.

    if exclude is not None:
        for i in exclude:
            num.remove(i)

    random.shuffle(num)
    # random.shuffle return None.
    # do not write `return random.shuffle(num)`

    return num
