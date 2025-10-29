"""Module for the Chain class and related helper functions."""

# pylint: disable=no-absolute-import
# Standard library imports
import re
import copy


# Local imports
from InDelsTopo.alphabet import Alphabet
from InDelsTopo.base import ChainBlockBase
from InDelsTopo.block import Block


# Auxiliary Functions
# ------------------------------------------------------------------------
def _combine_blocks_alphabet(blocks_list, alphabet=None):
    """Combine the alphabets of a collection of blocks into a single Alphabet.

    This function aggregates all unique symbols in the alphabets of the blocks in `W`,
    merges them with an optional initial alphabet, and returns a new
    `Alphabet` object containing the combined symbols.

    Args:
        W (iterable): A collection of Block objects.
        alphabet (Alphabet, optional): An initial `Alphabet` to merge with.
            If None, starts with an empty `Alphabet`.

    Returns:
        Alphabet: A new Alphabet instance containing all symbols
        from the provided blocks and the initial alphabet.
    """
    initial_alphabet = Alphabet([]) if alphabet is None else alphabet
    alphabets = [B.get_alphabet().letters for B in blocks_list]
    all_letters = initial_alphabet.letters
    for letters in alphabets:
        all_letters.update(letters)
    new_alphabet = Alphabet(list(all_letters.keys()))
    return new_alphabet


def _chain_str_constructor(expression, alphabet=None, prod_symbol=None):
    """Main function for constructing a chain from its string representation.

    Args:
        expression (str): The string representation of a chain. Each term can
            include an optional integer coefficient followed by a block, e.g.,
            '2x_0(1,a_1)x_1(1,a_2)x_2 - y_0(1,b_1)y_1'.
        alphabet (Alphabet, optional): An existing Alphabet object to use.
            Defaults to None, in which case a new empty Alphabet instance is
            created and populated with letters from the blocks.
        prod_symbol (str, optional): The product symbol to use in blocks.
            Must be one of {'', '*', '.'}. If None, it is inferred from the
            expression: '*' if '*' appears, '.' if '.' appears, or '' (concatenation)
            otherwise.

    Returns:
        tuple:
            list_coeffs (list of int): The coefficients of each term in the chain.
            blocks_x_factors_edges (list of Block): A list of Block objects
                corresponding to the parsed blocks. Each Block contains
                `x_factors` (list of int or SymPy expression) and `edges`
                (list of SymPy symbol).
            alphabet (Alphabet): The Alphabet instance used, either the provided
                one or a newly created and populated instance.
            prod_symbol (str): Either '', '.', or '*'.

    Notes:
        - The function standardizes the input string by converting '**' to '^'
          and removing whitespace.
        - Coefficients are parsed using `_convert_coeff`.
        - Blocks are parsed using the `Block` constructor, which handles
          internal x_factors and edges.
        - The product symbol is unified across all blocks for consistency.
    """

    alphabet = Alphabet([]) if alphabet is None else alphabet
    expression = expression.replace("**", r"^")
    if prod_symbol is None:
        if "*" in expression:
            prod_symbol = "*"
            prod_symbol_message = "*"
            # If '.' also included by mistake, replace with '*'
            expression = expression.replace(".", "*")
        elif "." in expression:
            prod_symbol = "."
            prod_symbol_message = "."
        else:
            prod_symbol = ""
            prod_symbol_message = "(concatenation)"
        print("Product symbol set to " + prod_symbol_message)

    expression = re.sub(r"\s+", "", expression)
    split_pattern = r"((?:^\+?\-?\d*)|(?:\+\d*)|(?:\-\d*))"
    terms_list = re.split(split_pattern, expression)
    terms_list.pop(0)
    list_coeffs = [_convert_coeff(c) for c in terms_list[::2]]
    list_blocks = terms_list[1::2]

    blocks_x_factors_edges = [
        Block(B, alphabet=alphabet, prod_symbol=prod_symbol) for B in list_blocks
    ]

    return list_coeffs, blocks_x_factors_edges, alphabet, prod_symbol


def _convert_coeff(string):
    """Convert a string representing a coefficient into an integer.

    Args:
        string (str): The input string representing a coefficient.

    Returns:
        int: The integer value of the coefficient.
    """
    if string in ["", "+", "-"]:
        string += "1"
    return int(string)


def _clean_dictionary_blocks(dict_blocks):
    """Remove empty or zero-valued blocks from a dictionary of blocks. This
    function is intended for internal use within the module.

    Args:
        dict_blocks (dict): A dictionary with Block objects as keys.

    Effects:
        Mutates `dict_blocks` in place by deleting:
          - Any keys whose value equals 0.
          - All occurrences of Block(0).
    """
    list_blocks = copy.deepcopy(list(dict_blocks.keys()))
    for block in list_blocks:
        if dict_blocks[block] == 0:
            del dict_blocks[block]
    # Delete '0' block
    while Block(0) in dict_blocks:
        del dict_blocks[Block(0)]


# ------------------------------------------------------------------------


# Chain Class
class Chain(ChainBlockBase):
    """
    Represents a chain in the Insertion Chain Complex:
    a linear combination of valid blocks with integer coefficients.

    It can be initialized from a string expression, a list of coefficients and blocks,
    or a dictionary mapping blocks to coefficients.

    Supports algebraic operations (addition, subtraction, and scalar multiplication),
    equality checking, string and LaTeX representations, and computation of the boundary.

    Attributes:
        dim (int): Maximum dimension among the blocks in the chain, or -1 if empty.

    Internal_attributes:
        _expression(SymPy expression): represents the chain as a sum of block expressions.
        _alphabet (Alphabet): The Alphabet instance used for all blocks in the chain.
        _dict_blocks (dict): (for internal use only)
            dictionary mapping Block objects to integer coefficients.
        _prod_symbol (str): Either '', '*', or '.'.
    """

    def __init__(
        self,
        expression=None,
        alphabet=None,
        prod_symbol=None,
        *,
        list_coeffs=None,
        list_blocks=None,
        dict_blocks=None,
    ):
        """
        Represents a chain in the Insertion Chain Complex:
        a linear combination of valid blocks with integer coefficients.

        Each term can include an integer coefficient followed by a block, e.g.,
        '2abc(1,a)b^ac(1,b)ac - (1,c)(1,a)b^2'.

        The class can be initialized from:
            - a string expression (preferred way for end users),
            - a list of coefficients and blocks, or
            - a dictionary mapping blocks to coefficients.

        Args:
            expression (str, optional): Chain written as a string of blocks. This is the
                preferred way for end users to construct a chain.
            alphabet (Alphabet, optional): Alphabet object defining the set of letters.
                If set to None, a new empty alphabet is created.
            prod_symbol (str, optional): Product symbol to use in blocks.
                Must be one of {'', '*', '.'}. If None, it is inferred from the expression:
                '*' if '*' appears, '.' if '.' appears, or '' otherwise.
            list_coeffs (list of int, optional): Integer coefficients of the chain terms.
                Must be provided together with `list_blocks`.
                Ignored if `expression` or `dict_blocks` is given.
            list_blocks (list of Block, optional): Block objects corresponding to the coefficients
                in `list_coeffs`. Ignored if `expression` or `dict_blocks` is given.
            dict_blocks (dict, optional): Dictionary mapping Block objects to integer coefficients.
                If provided, it takes precedence over `expression`,`list_coeffs`, and `list_blocks`.
        """
        super().__init__()
        alphabet = Alphabet([]) if alphabet is None else alphabet

        # If dict_blocks provided use that one
        if dict_blocks is not None:
            # Check all keys are instances of Block
            if not all(isinstance(b, Block) for b in dict_blocks.keys()):
                raise TypeError("All keys of dict_blocks must be instances of Block")

            # Check all values are integers
            if not all(isinstance(c, int) for c in dict_blocks.values()):
                raise TypeError("All values of dict_blocks must be integers")

            self._dict_blocks = dict_blocks
            self._alphabet = _combine_blocks_alphabet(dict_blocks.keys(), alphabet)

            # Use provided prod_symbol or use a suitable: concatenation for symbols of
            # length 1, "." otherwise
            if prod_symbol in ["", ".", "*"]:
                self._prod_symbol = prod_symbol
            elif (
                max([len(symbol) for symbol in self._alphabet.letters_str], default=1)
                > 1
            ):
                self._prod_symbol = "."
            else:
                self._prod_symbol = ""

        # Otherwise, if expression is present, or list_coeff, list_blocks are
        # present, use that.
        else:
            if isinstance(expression, str):
                list_coeffs, list_blocks, alphabet, prod_symbol = (
                    _chain_str_constructor(expression, alphabet, prod_symbol)
                )
                self._alphabet = alphabet
                self._prod_symbol = prod_symbol

            else:  # Check list_coeffs and list_blocks are valid and create alphabet and prod_symbol

                # Check list_coeffs and list_blocks have same length
                if len(list_coeffs) != len(list_blocks):
                    raise ValueError(
                        "list_coeffs and list_blocks must have the same length"
                    )

                # Check all elements in list_blocks are instances of Block
                if not all(isinstance(b, Block) for b in list_blocks):
                    raise TypeError(
                        "All elements of list_blocks must be instances of Block"
                    )

                # Check all elements in list_coeffs are integers
                if not all(isinstance(c, int) for c in list_coeffs):
                    raise TypeError("All elements of list_coeffs must be integers")

                # Set alphabet

                self._alphabet = _combine_blocks_alphabet(dict_blocks.keys(), alphabet)
                # Use provided prod_symbol or use a suitable: concatenation for symbols of
                # length 1, "." otherwise
                if prod_symbol in ["", ".", "*"]:
                    self._prod_symbol = prod_symbol
                elif (
                    max(
                        [len(symbol) for symbol in self._alphabet.letters_str],
                        default=1,
                    )
                    > 1
                ):
                    self._prod_symbol = "."
                else:
                    self._prod_symbol = ""

            # Create dictionary
            self._dict_blocks = {}
            for i, block in enumerate(list_blocks):
                if block not in self._dict_blocks:
                    self._dict_blocks[block] = list_coeffs[i]
                else:
                    self._dict_blocks[block] += list_coeffs[i]

        # Clean dict_blocks
        _clean_dictionary_blocks(self._dict_blocks)

        # Set up attributes
        self._create_expression()
        self.dim = max([B.dim for B in self._dict_blocks] + [-1])

    def _create_expression(self):
        """Build the SymPy expression from the chain's blocks and
        coefficients."""
        # Create expression
        self._expression = int(0) + sum(
            [
                int(self._dict_blocks[block]) * block.get_expression()
                for block in self._dict_blocks
            ]
        )

    def __add__(self, other):
        new_dict = copy.deepcopy(self._to_chain()._dict_blocks)
        other = other._to_chain()
        for block in other._dict_blocks:
            if block in new_dict:
                new_dict[block] += other._dict_blocks[block]
            else:
                new_dict[block] = other._dict_blocks[block]
        _clean_dictionary_blocks(new_dict)
        return Chain(dict_blocks=new_dict)

    def __rmul__(self, coeff):
        dict_blocks = self._to_chain()._dict_blocks.copy()
        if coeff != 0:
            new_dict = {c: coeff * dict_blocks[c] for c in dict_blocks}
        else:
            new_dict = {}
        return Chain(
            dict_blocks=new_dict, alphabet=self._alphabet, prod_symbol=self._prod_symbol
        )

    def boundary(self):
        """Return the boundary of the chain as a new Chain."""
        new_dict = {}
        for block in self._dict_blocks:
            coeff = self._dict_blocks[block]
            dict_boundary = block.boundary().get_dict_blocks()
            for new_block in dict_boundary:
                if new_block in new_dict:
                    new_dict[new_block] += coeff * dict_boundary[new_block]
                else:
                    new_dict[new_block] = coeff * dict_boundary[new_block]

        _clean_dictionary_blocks(new_dict)
        return Chain(
            dict_blocks=new_dict, alphabet=self._alphabet, prod_symbol=self._prod_symbol
        )

    def get_dict_blocks(self):
        """Return a copy of the internal blocks dictionary."""
        return self._dict_blocks.copy()
