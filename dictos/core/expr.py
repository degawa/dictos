import sympy as sp
from sympy.printing.str import StrPrinter
from sympy.printing.latex import LatexPrinter
from sympy.printing.repr import ReprPrinter


class Expr(sp.Expr):
    """
    A class extending SymPy's Expr to preserve expression display order.

    This class wraps SymPy expressions while maintaining the original order
    of terms in the expression, which is useful for cases where the display
    order matters (e.g., educational purposes or specific notation requirements).

    Attributes:
        arg (sympy.Expr): The underlying SymPy expression

    Example:
        >>> import dictos
        >>> stencil = [-2, -1, 0, 1, 2]
        >>> c = dictos.finite_difference.equation(stencil)
        >>> str(c)
        '(-8*f_{-1} + f_{-2} + 8*f_{1} - f_{2})/(12*h)' # Sorted by sympy while printing
        >>> d = dictos.Expr(c)
        >>> str(d)
        '(f_{-2} - 8*f_{-1} + 8*f_{1} - f_{2})/(12*h)' # Original order is preserved
    """

    def __new__(cls, arg):
        """
        Create a new Expr instance.

        This method ensures proper initialization of the SymPy expression
        and handles cases where the input is already an Expr instance.

        Args:
            arg (sympy.Expr): A SymPy expression to be wrapped

        Returns:
            Expr: A new Expr instance wrapping the input expression

        Example:
            >>> import dictos
            >>> from sympy.abc import x
            >>> expr = dictos.Expr(x+1)
            >>> isinstance(expr, dictos.Expr)
            True
        """

        if isinstance(arg, Expr):
            return arg
        obj = super().__new__(cls, arg)
        obj.arg = arg
        return obj

    def __str__(self) -> str:
        """
        Return a string representation of the Expr.

        This method uses StrPrinter with 'none' ordering to convert
        the Expr to conform to existing notaiton.

        Returns:
            str: A human-readable string representation

        Example:
            >>> import dictos
            >>> from sympy.abc import x, y
            >>> f = dictos.Expr(x + 2*y)
            >>> str(f)
            'x + 2*y'
        """
        return StrPrinter({"order": "none"}).doprint(self.arg)

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the Expr.

        This method provides an unambiguous representation of the expression
        that can be used to recreate the object. It uses ReprPrinter with
        'none' ordering and wraps the result in Expression().

        Returns:
            str: A string that can be evaluated to recreate the object

        Example:
            >>> import dictos
            >>> from sympy.abc import x, y
            >>> f = dictos.Expr(x + 2*y)
            >>> repr(f)
            "Expression(Add(Symbol('x'), Mul(Integer(2), Symbol('y'))))"
        """
        sympy_repr = ReprPrinter({"order": "none"}).doprint(self.arg)
        return f"Expression({sympy_repr})"

    def _latex(self) -> str:
        """Generate a LaTeX representation of the Expr.

        This method creates a LaTeX string using LatexPrinter with 'none' ordering,
        preserving the original term order. This method is primarily used
        by Jupyter notebooks and other LaTeX-aware display systems.

        Returns:
            str: LaTeX representation of the expression

        Note:
            This is an internal method typically called by LaTeX rendering
            systems. Direct use is discouraged - use display() or IPython's
            display mechanisms instead.

        Example:
            >>> import dictos
            >>> from sympy.abc import x, y
            >>> f = dictos.Expr(x + 2*y)
            >>> f._latex()
            'x + 2y'
        """
        return LatexPrinter({"order": "none"}).doprint(self.arg)

    def toSympyExpr(self) -> sp.Expr:
        """Convert to a standard SymPy expression.

        This method extracts the underlying SymPy expression, which can be useful when
        you need to use SymPy's standard functionality that expects a
        regular SymPy Expr instance.

        Returns:
            sympy.Expr: The underlying SymPy expression

        Example:
            >>> import dictos
            >>> import sympy as sp
            >>> from sympy.abc import x
            >>> expr = dictos.Expr(x + 1)
            >>> sympy_expr = expr.toSympyExpr()
            >>> isinstance(sympy_expr, sp.Expr)
            True
        """
        return self.arg
