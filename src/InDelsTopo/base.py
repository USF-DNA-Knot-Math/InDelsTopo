"""Module for base class for Chain and Block."""

# Third party imports
import sympy as sym


# Chain Class
class ChainBlockBase:
    """Base class representing a Chain or a Block in the InDelsTopo package.

    This class provides the core functionality for:
        - algebraic operations on Chains and Blocks (addition, subtraction, scalar multiplication)
        - representations (__str__, __repr__, LaTeX)
        - access to underlying expression and Alphabet
        - retrieval of constituent blocks

    Attributes:
        _alphabet (Alphabet or None): The underlying alphabet of the Chain/Block.
        _prod_symbol (str or None): Symbol used for multiplication in the expression.
        _expression (SymPy expression or None): Symbolic expression representing the Chain/Block.
        dim (int or None): dimension of the Chain/Block.
        _dict_blocks: Dictionary of blocks.
    """

    def __init__(self, alphabet=None, prod_symbol=None, expression=None, dim=None):
        self._alphabet = alphabet
        self._prod_symbol = prod_symbol
        self._expression = expression
        self.dim = dim
        self._dict_blocks = {}

    def _to_chain(self):
        """Returns self.

        Useful for adding/substracting Blocks. Intended to be overridden
        by Block.
        """
        return self

    def _repr_latex_(self):
        """Return the LaTeX representation of the Chain or Block for Jupyter
        display."""
        dict_symbols = {"": "", ".": "\\cdot ", "*": "*"}
        try:
            return (
                "$"
                + sym.latex(
                    self._expression, mul_symbol=dict_symbols[self._prod_symbol]
                )
                + "$"
            )
        except Exception:  # pylint: disable=W0703
            return str(self._expression)

    def __str__(self):
        try:
            return (
                str(self._expression)
                .replace("**", r"^")
                .replace("*", self._prod_symbol)
            )
        except Exception:  # pylint: disable=W0703
            return str(self._expression)

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        if not hasattr(other, "_to_chain"):
            return NotImplemented
        return self._to_chain() + other._to_chain()

    def __rmul__(self, coeff):
        if not isinstance(coeff, int):
            return NotImplemented
        return coeff * self._to_chain()

    def __sub__(self, other):
        if not hasattr(other, "_to_chain"):
            return NotImplemented
        return self._to_chain() + (-1) * other._to_chain()

    def __neg__(self):
        return (-1) * self._to_chain()

    def __eq__(self, other):
        if not hasattr(other, "_expression"):
            return NotImplemented
        return self._expression == other._expression

    def get_blocks(self):
        """Return a list of the blocks in the chain."""
        return list(self._to_chain().get_dict_blocks().keys())

    def get_expression(self):
        """Return a SymPy expression representing this Chain."""
        return self._expression

    def get_alphabet(self):
        """Returns the underlying Alphabet object."""
        return self._alphabet

    def __hash__(self):
        return self._expression.__hash__()

    def get_dict_blocks(self):
        """Returns the block dictionary"""
        return self._to_chain().get_dict_blocks()
