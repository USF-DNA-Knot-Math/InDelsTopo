"""Module for the Chain class and related helper functions."""

# pylint: disable=no-absolute-import
# Standard library imports
import re
import warnings

# Import external packages
import sympy as sym


# Local imports
from InDelsTopo.alphabet import Alphabet
from InDelsTopo.base import ChainBlockBase
from InDelsTopo.utils import _expand_symp_word, _powerset

# Auxiliary Functions
# -----------------------------------------------------------------------------------------


def _is_subsequence(subseq, seq):
    """Check whether `subseq` is a subsequence of `seq`."""
    # Adapted from https://stackoverflow.com/a/24017747
    iterator = iter(seq)
    return all(x in iterator for x in subseq)


def _block_str_constructor(expression, alphabet=None, prod_symbol=None):
    """Main function for creating a block, given its representation as
    "x_0(1,a_1)x_1(1,a_2)...(1,a_m)x_m".

    Handles powers using '^' or '**' (e.g., 'a^2', 'b**3') and product symbols
    ('*', '.', or concatenation).

    Args:
        expression (str): The string representation of a block.
        alphabet (Alphabet, optional): An existing Alphabet object to use.
            Defaults to None, in which case a new empty Alphabet instance is
            created and populated with letters from the expression.
        prod_symbol (str, optional): The product symbol to use. Must be one of
            {'', '*', '.'}. If None, it is inferred from the expression:
            set to '*' if '*' appears, to '.' if '.' appears, or to '' (concatenation)
            otherwise.

    Returns:
        tuple:
            x_factors (list of SymPy Expressions or int): A list of factors
                corresponding to the x_i components (products of alphabet letters
                with their powers).
            edges (list of SymPy symbols): A list of edge elements
                (letters from the alphabet) extracted from the (1,a_i) terms of the block.
            alphabet (Alphabet): The Alphabet instance used, either the provided
                one or a newly created and populated instance.
            prod_symbol (str): Either '', '.', or '*'.

    Raises:
        ValueError: If an edge in the expression has an exponent greater than 1,
            indicating an invalid block.

    Notes:
        This function standardizes the input string by converting '**' to '^',
        unifying product symbols, and parsing alternating sequences of x-factors
        and edges. It is tolerant of mixed product symbols but will enforce
        consistency for subsequent parsing.
    """
    alphabet = Alphabet([]) if alphabet is None else alphabet
    # Make sure all powers and products are given by the same symbol,
    # Especially if the output of a str() from SymPy has been applied.
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

    # Separate expression into the x_factors and edges
    pattern = r"(\(1,[a-zA-Z0-9_]+\)(?:\^\d+)?|[a-zA-Z0-9_\^" + prod_symbol + "]+)"
    list_factors = re.split(pattern, expression)
    list_factors = [factor for factor in list_factors if factor != ""]

    # x_factors
    if prod_symbol == "":
        pattern_alpha = r"(?P<base>[a-zA-Z0-9_])(?:\^(?P<power>\d+))?"
    else:
        pattern_alpha = r"(?P<base>[a-zA-Z0-9_]+)(?:\^(?P<power>\d+))?"
    re_alpha = re.compile(pattern_alpha)

    pattern_edge = r"\(1,(?P<base>[a-zA-Z0-9_]+?)\)(?:\^(?P<power>\d+))?"
    re_edge = re.compile(pattern_edge)

    # Find sequences of edges and x_factors
    currently_at_alpha = True
    current_factor = 1

    x_factors = []
    edges = []
    for factor in list_factors:
        if currently_at_alpha:
            # Check if we have an edge or not
            match_edge = re_edge.match(factor)
            # If we hit an edge, the current factor is the previous alpha.
            if match_edge is not None:
                x_factors.append(current_factor)
                current_factor = 1
                currently_at_alpha = False
            else:  # We got something to add to alpha
                new_factors = re_alpha.findall(factor)
                for base, power in new_factors:
                    alphabet.add_letter(base)
                    if power == "":
                        power = 1
                    else:
                        power = int(power)
                    current_factor *= alphabet.letters[base] ** power
        if not currently_at_alpha:
            edge = re_edge.match(factor)
            if not (edge["power"] is None) and (int(edge["power"]) > 1):
                # If the power of an edge is not 1, we have an invalid block, thus we get
                # the zero block.
                return [0], [], alphabet, prod_symbol
            alphabet.add_letter(edge["base"])
            edges.append(alphabet.letters[edge["base"]])
            currently_at_alpha = True
    x_factors.append(current_factor)
    return x_factors, edges, alphabet, prod_symbol


# -----------------------------------------------------------------------------------------
# Block Class


class Block(ChainBlockBase):
    """Represents a Block in the Insertion Chain Complex.

    A Block can be a valid block or "0" if the provided data corresponds to an invalid block.
    It has the formal expression x_0(1,a_1)x_1...(1,a_m)x_m, where the a_i are single symbols
    from an alphabet, and the x_i are words over that alphabet (including possibly the empty word).
    Internally, blocks are stored as a list of factors x_0, ..., x_k and a list of edge
    symbols a_1, ..., a_k, represented as SymPy expressions and symbols.

    Blocks can be initialized from a string expression or directly via x_factors and edges.

    Attributes:
        dim (int): Number of edges (dimension) of the block.
        max_word (SymPy expression): Maximal word of the block.
        min_word (SymPy expression): Minimal word of the block.
        # _x_factors (list of SymPy expressions or int): Factors of the block.
        # _edges (list of SymPy symbols): Symbols forming the edges of the block.
        # _alphabet (Alphabet): Alphabet used for letters in the block.
        # _expression (SymPy expression): SymPy product representing the block.
    """

    def __init__(
        self,
        expression=None,
        alphabet=None,
        prod_symbol=None,
        *,
        x_factors=[int(0)],
        edges=None,
    ):
        """Initializes a Block.

        A Block can be either a valid block or the zero block. It is internally
        represented by a list of factors `_x_factors` and a list of edge symbols
        `_edges`.

        The block can be initialized from:
            - A string `expression` in the form "x_0(1,a_1)x_1...(1,a_m)x_m"
              (handling powers '^' or '**' and product symbols '*', '.', or concatenation),
            - Directly from `_x_factors` and `_edges` lists,

        Args:
            expression (str, optional): String representation of the block.
            alphabet (Alphabet, optional): Alphabet defining the letters.
                Defaults to an empty Alphabet.
            prod_symbol (str, optional): Product symbol to use ('', '*', or '.').
                Inferred from the expression if None.
            x_factors (list of SymPy expressions or int, optional):
                Factors x_0, ..., x_k of the block.
            edges (list of SymPy symbols, optional): Edge symbols a_1, ..., a_k
                corresponding to the block.
        """
        super().__init__()  # initialize base class attributes
        x_factors = [int(0)] if x_factors is None else x_factors
        edges = [] if edges is None else edges
        alphabet = Alphabet([]) if alphabet is None else alphabet

        # Use expression if provided
        if isinstance(expression, str):
            if expression == "1":
                x_factors = [int(1)]
                edges = []
                prod_symbol = ""
            else:
                x_factors, edges, alphabet, prod_symbol = _block_str_constructor(
                    expression, alphabet, prod_symbol
                )

        if len(x_factors) == 0:
            self._x_factors = [int(0)]
            self._edges = []
            self.dim = -1

        elif len(x_factors) - len(edges) != 1:
            raise Exception("x_factors must be one unit longer than edges.")

        # Set attributes
        self._x_factors = x_factors
        self._edges = edges
        self.dim = len(edges)
        self._prod_symbol = prod_symbol
        self._alphabet = alphabet

        # Get canonical form
        self._canonical_form()

        # Check is a valid block
        if not self._is_valid_block():
            self._x_factors = [int(0)]
            self._edges = []
            self.dim = -1

        # Expression
        edges = [self._alphabet.get(str(a), 1) for a in self._edges] + [sym.sympify(1)]
        self._expression = sym.prod(
            [val for pair in zip(self._x_factors, edges) for val in pair]
        )

        # Maximal and minimal words
        if self.dim < 0:
            self.max_word = sym.sympify(0)
            self.min_word = sym.sympify(0)
        elif self.dim == 0:
            self.max_word = self._expression
            self.min_word = self._expression
        else:
            self.max_word = self.get_vertex(list(range(1, self.dim + 1)), as_word=True)
            self.min_word = self.get_vertex([], as_word=True)

    def _canonical_form(self):
        """Gets the canonical form of the block by transforming each factor of
        the form '(1,a)a^r' into 'a^r(1,a)'."""
        for k in range(self.dim, 0, -1):
            factors = list(sym.sympify(self._x_factors[k]).as_coeff_mul()[1])
            try:
                first_factor = factors.pop(0)
            except IndexError:
                first_factor = sym.sympify(1)
            if first_factor.as_base_exp()[0] == self._edges[k - 1]:
                self._x_factors[k - 1] = self._x_factors[k - 1] * first_factor
                self._x_factors[k] = sym.Mul(sym.prod(factors))

    def _is_valid_block(self):
        """Check whether a block, already in canonical form, is valid."""
        if self.dim > 0:
            position_ones = [
                i for i, s in enumerate(self._x_factors) if s == sym.sympify(1)
            ]
            for k in position_ones:
                if 0 < k < self.dim and self._edges[k] == self._edges[k - 1]:
                    return False
        elif self.dim == 0:
            if self._x_factors[0] == int(0):
                return False
        return True

    def _to_chain(self):
        """Convert the block to a Chain object.

        Useful to enable algebraic operations with chains.
        """
        from InDelsTopo.chain import Chain  # pylint: disable=import-outside-toplevel

        if self.dim >= 0:
            return Chain(dict_blocks={self: 1})
        return Chain(dict_blocks={})

    def get_face(self, indices_plus, indices_minus):
        """Return the face σ(indices_plus, indices_minus) of the block σ.

        indices_plus and indices_minus are disjoint subsets of {1, ..., m}, where m is the dimension
        of the block. They must be given as lists of integers, and may be empty.
        The output is the block σ(indices_plus, indices_minus), which is either a valid
        Block of dimension m − |indices_plus ∪ indices_minus| or the zero Block.

        Args:
            indices_plus (list of int): Indices in {1, ..., m} of edges (1, a_i) collapsed to a_i.
            indices_minus (list of int): Indices in {1,..., m} of edges (1, a_i) collapsed to 1.

        Returns:
            face (Block): The resulting face block σ(indices_plus, indices_minus) if valid, 
                or the zero Block, otherwise.
        """
        indices_all = indices_plus + indices_minus

        if len(set(indices_plus).intersection(indices_minus)) > 0:
            raise ValueError("indices_plus and indices_minus must be disjoint")

        if len(set(indices_all).difference(set(range(1, self.dim + 1)))) > 0:
            warnings.warn(
                "indices_plus or indices_minus contains elements outside of {1, ..."
                + str(self.dim)
                + "} which have been ignored"
            )

        new_x_factors = [self._x_factors[0]]
        new_edges = []
        for i in range(1, self.dim + 1):
            if i in indices_all:
                current_x_factor = new_x_factors[-1]
                if i in indices_plus:
                    x = self._edges[i - 1]
                else:
                    x = sym.sympify(1)
                new_x_factors[-1] = current_x_factor * x * self._x_factors[i]
            else:
                new_edges.append(self._edges[i - 1])
                new_x_factors.append(self._x_factors[i])
        return Block(
            x_factors=new_x_factors,
            edges=new_edges,
            alphabet=self._alphabet,
            prod_symbol=self._prod_symbol,
        )

    def get_upper_facets(self):
        """Return the list of upper facets of the block.

        Each upper facet is obtained by collapsing a single edge (1, a_i) to a_i,
        for i = 1, ..., m, where m is the dimension of the block.

        Returns:
            upper_facets (list[Block]): A list of Blocks representing the upper facets of the block.
        """
        upper_facets = [self.get_face([i], []) for i in range(1, self.dim + 1)]
        return upper_facets

    def get_lower_facets(self):
        """Return the list of lower facets of the block.

        Each lower facet is obtained by collapsing a single edge (1, a_i) to 1,
        for i = 1, ..., m, where m is the dimension of the block. Invalid blocks
        are excluded from the output.

        Returns:
            list of Block: A list of valid Blocks representing the lower facets of the block.
        """
        lower_facets = [self.get_face([], [i]) for i in range(1, self.dim + 1)]
        lower_facets = [C for C in lower_facets if C.dim >= 0]
        return lower_facets

    def get_all_facets(self):
        """Return all facets of the block. Combines both the upper and lower
        facets of the block.

        Returns:
            list of Block: A list of Blocks representing all facets of the block.
        """
        return self.get_upper_facets() + self.get_lower_facets()

    def get_all_faces(self, include_self=False):
        """Return all faces of the block, ordered by dimension.

        Faces are obtained recursively by taking all facets of the block,
        then all facets of those facets, and so on, down to dimension 0.

        Args:
            include_self (bool, optional): If True, include the block itself
                in the returned list. Defaults to False.

        Returns:
            list of Block: A list of Blocks representing all faces, orted by increasing dimension.
        """
        if include_self:
            all_faces = set([self])
        else:
            all_faces = set()
        current_faces = [self]
        for _ in range(self.dim):
            current_faces = [
                facet for c in current_faces for facet in c.get_all_facets()
            ]
            all_faces = all_faces.union(current_faces)
        all_faces = list(all_faces)
        all_faces.sort(key=lambda x: x.dim)
        return all_faces

    def get_vertex(self, indices, as_word=True):
        """Return the vertex v_I(σ) of the block determined by a sequence of
        I = indices.

        The vertex is obtained by collapsing the edges (1,a_i) indexed by 'I' into a_i,
        and the remaining ones to 1. Returns the vertex as a 0-block,
        or as a SymPy expression if as_word=True.

        Args:
            indices (list of int): Indices in {1,..., dim(σ)}
            as_word (bool, optional): If True, return the vertex as a SymPy
                expression (word). If False, return the corresponding Block.
                Defaults to True.

        Returns:
            Block or SymPy Expression: The resulting vertex as a Block or Sympy expression,
            depending on `as_word`.
        """

        if len(set(indices).difference(set(range(1, self.dim + 1)))) > 0:
            warnings.warn(
                "indices contains elements outside of {1, ..."
                + str(self.dim)
                + "} which have been ignored"
            )

        indices_minus = list(set(range(1, self.dim + 1)).difference(indices))
        vertex = self.get_face(indices, indices_minus)
        if as_word:
            return vertex.get_expression()
        return vertex

    def get_all_vertices(self, as_words=True):
        """Return all vertices of the block.

        Args:
            as_words (bool, optional): If True, return vertices as SymPy expressions
                (words). If False, return them as Block objects. Defaults to True.

        Returns:
            list of Block or SymPy Expression: All vertices of the block as Blocks
            or words, depending on `as_words`.
        """
        vertices = set()
        all_indices = _powerset(list(range(1, self.dim + 1)))
        for indices in all_indices:
            vertices.add(self.get_vertex(list(indices), as_words))
        return list(vertices)

    def boundary(self):
        """Computes the boundary of a block."""
        from InDelsTopo.chain import Chain  # pylint: disable=import-outside-toplevel

        result = Chain(dict_blocks={})
        for i in range(1, self.dim + 1):
            result += (-1) ** (i + 1) * (
                self.get_face([i], []) - self.get_face([], [i])
            )
        return result

    def _le_recursive(self, other_block):
        """Recursively determine whether this block is a face of another block.

        This internal helper performs the recursive part of the face relation
        once trivial and extremal cases have been handled. It is not meant to be
        called directly—use the operator form `A <= B` instead.

        Args:
            other_block (Block): The block to compare against.

        Returns:
            bool: True if `self` is a (possibly proper) face of `other_block`,
            False otherwise.
        """
        # If same dimension, compare block equality
        if self.dim == other_block.dim:
            return self == other_block

        # Compare edges sequences pylint: disable=protected-access
        if not _is_subsequence(self._edges, other_block._edges):
            return False

        # Expand and compare the maximal words symbol by symbol
        word1 = self.max_word
        word2 = other_block.max_word
        word1_extended = [
            symbol
            for symbol, times in _expand_symp_word(word1)[0]
            for i in range(times)
        ]
        word2_extended = [
            symbol
            for symbol, times in _expand_symp_word(word2)[0]
            for i in range(times)
        ]
        if not _is_subsequence(word1_extended, word2_extended):
            return False

        # Recursive step
        # If maximal words are equal, descend through upper facets;
        # otherwise, through lower facets.
        if word1 == word2:
            return any(
                self._le_recursive(facet) for facet in other_block.get_upper_facets()
            )
        else:
            return any(
                self._le_recursive(facet) for facet in other_block.get_lower_facets()
            )

    def __le__(self, other_block):
        """
        Determine whether this Block is a face (≤) of another Block.

        A Block `A` is considered less than or equal to another Block `B` if it
        represents a lower-dimensional face of `B`.
        """
        # Make sure both are block instances
        if not isinstance(other_block, Block):
            return NotImplemented

        # "Empty" block is face of all blocks
        if self.dim == -1:
            return True

        # Compare dimensions
        if self.dim > other_block.dim:
            return False

        # Apply the comparison recursively
        return self._le_recursive(other_block)

    def __lt__(self, other_block):
        """Return True if this Block is a proper face (<) of another Block."""
        if not isinstance(other_block, Block):
            return NotImplemented
        if self.__eq__(other_block):
            return False
        return self.__le__(other_block)
