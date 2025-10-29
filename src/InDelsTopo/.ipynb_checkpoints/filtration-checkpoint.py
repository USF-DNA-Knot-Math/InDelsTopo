"""Module for the Filtrationclass and related helper functions."""

# pylint: disable=no-absolute-import
# pylint: disable=consider-using-from-import
# Standard library imports
import re

# import warnings

# Import external packages
import sympy as sym
import numpy as np


# Local imports
from InDelsTopo.block import Block
from InDelsTopo.alphabet import Alphabet
from InDelsTopo.utils import _combine_blocks_alphabet, _expand_symp_word
import InDelsTopo.graphics as graphics


# Auxiliary functions
# -----------------------------------------------------------------------------------------
def _compute_new_index_lower_facet(j, i, exp_i, base_before, base_after):
    """
    Compute the updated index for a lower facet of a block.

    Determines how the index `j` changes when removing the `i`-th edge
    from a block, depending on the exponent and the surrounding symbols.

    Args:
        j (int): The original index.
        i (int): The index of the edge being removed.
        exp_i (int): The exponent of the symbol at position `i`.
        base_before (SymPy expression): Symbol before position `i`.
        base_after (SymPy expression): Symbol after position `i`.

    Returns:
        int: The updated index corresponding to `j`.
    """
    if j < i:
        return j
    if exp_i > 1:
        return j
    if base_before != base_after:
        return j - 1
    return j - 2


def _convert_words_to_blocks(words, prod_symbol=None, alphabet=None):
    """
    Convert a list of strings into Block objects.
    Args:
        words (list of str): The collection of words to convert.
        prod_symbol (str, optional): The symbol used for product/concatenation.
            If None, the symbol is inferred from the input words.
        alphabet (Alphabet, optional): Alphabet object defining valid symbols.
            If None, a new empty `Alphabet` is created.

    Returns:
        list of str: The normalized list of words, ready for block conversion.
    """
    alphabet = Alphabet([]) if alphabet is None else alphabet
    words = [w.replace("**", r"^") for w in words]
    all_combined = " ".join(words)
    if prod_symbol is None or not (prod_symbol in ["", ".", "*"]):
        if "*" in all_combined:
            prod_symbol = "*"
            prod_symbol_message = "*"
            # If '.' also included by mistake, replace with '*'
            words = [w.replace(".", "*") for w in words]
        elif "." in all_combined:
            prod_symbol = "."
            prod_symbol_message = "."
        else:
            prod_symbol = ""
            prod_symbol_message = "(concatenation)"
        print("Product symbol set to " + prod_symbol_message)

    words = [re.sub("\\s+", "", w) for w in words]
    blocks = [Block(w, alphabet=alphabet, prod_symbol=prod_symbol)
              for w in words]

    return blocks, alphabet, prod_symbol


def _candidate_indices_next(indices_lists):
    """
    Internal helper function used exclusively by `compute_d_skeleton` to identify
    all (k+1)-index sets whose k-index subsets appear in `indices_lists`.

    This function reduces the number of candidate blocks with a given maximal block
    that need to be checked for support in the vertex set.

    Args:
        indices_lists (list of list of int): List of index lists, each corresponding to
        k indices [i_1, ..., i_k] such that the maximal word with each indices set is a valid
        block in the filtration.

    Returns:
        list of list of int: List of new index lists representing (k+1)-indices.
    """

    k = len(indices_lists[0])  # length of indices sets.
    # Sort indices and make sure we dont have multiple copies
    indices_lists = list({tuple(sorted(indices)) for indices in indices_lists})

    # Make them into sets
    indices_lists = [set(indices) for indices in indices_lists]

    # Finding intersections
    n = len(indices_lists)
    intersections = [{} for i in range(n)]
    for i in range(n - 1):
        indices1 = indices_lists[i]
        for j in range(i + 1, n):
            indices2 = indices_lists[j]
            if len(indices1.intersection(indices2)) == k - 1:
                new_index = indices2.difference(indices1).pop()
                if new_index in intersections[i]:
                    intersections[i][new_index] += 1
                else:
                    intersections[i][new_index] = 1

    new_indices_list = []
    # Found Groups
    for i in range(n):
        counts = intersections[i]
        for new_index in counts:
            if counts[new_index] == k:
                indices = indices_lists[i].copy()
                indices.add(new_index)
                new_indices_list.append(sorted(list(indices)))
    return new_indices_list


# ------------------------------------------------------------------------------------------------


# Filtration Class
class Filtration:
    """
    Represents a lower-star filtration on a Letter Insertion Chain Complex C[W],
    for a set of words W with associated heights.

    It stores the blocks of the complex at each dimension with their corresponding heights.
    Provides methods for computing the Euler Characteristic Curve,
    persistent homology barcodes (with Z/Z2 coefficients),
    and a graphical representation of (accurate for low dimensions).

    One can access the k-th dimensional blocks by indexing the filtration as 'K[k]'.
    This returns a dictionary mapping each k-block to its height.

    Attributes:
        dim (int): Maximum dimension of the filtration.
        filtration_dict (dict): Maps dimension d to a dictionary of blocks and their heights.
        filtration_values (list): Sorted list of heights used in the filtration.

    Internal Attributes:
        _alphabet (Alphabet): Alphabet containing all symbols in W.
        _prod_symbol (str): Product symbol used for constructing blocks
            ('*', '.', or '' for concatenation).
        _positions_dict (dict or None): Stores vertex positions for graph representations.

    Notes:
        - It initializes as an empty filtration. It can be made into the filtration of insertion chain
            complexes of a set of words `W` by using the method `compute_d_skeleton(W, heights)`.
        - Blocks can be added or removed by using the methods `add_blocks` and `remove_blocks`.

    Example:
    >>> W = ["ab", "aab", "abb"]
    >>> heights = [0.1, 0.3, 0.5]
    >>> K = Filtration() # Creates an empty complex
    >>> K.compute_d_skeleton(W, heights, max_dim=5) # makes K = a filtration of C[W]
    >>> K[1]  # Access 1-dimensional blocks and their heights
    {a(1,a)b: 0.3, ab(1,b): 0.5}
    """

    def __init__(self, alphabet=None, prod_symbol=None):
        """
        Returns an empty Filtration object.

        Args:
            alphabet (Alphabet, optional): Alphabet object defining valid symbols.
                If None, a new empty `Alphabet` is created.
            prod_symbol (str, optional): The product symbol to use. Must be one of
                {'', '*', '.'}. If None, it is inferred from the expression:
                set to '*' if '*' appears, to '.' if '.' appears, or to '' (concatenation)
                otherwise.
        """
        # User inputs
        self._alphabet = alphabet
        self._prod_symbol = prod_symbol

        # Filtration data
        self.dim = -1
        self.filtration_dict = {}
        self.filtration_values = []

        # Used when creating a graph
        self._positions_dict = None

    def _facets_maximal_word(self, word, indices, which="all", casted=False):
        """
        Generates the (upper/lower/all) facets of a block given as maximal word-indices pair.

        Args:
            w (str or SymPy expression): The maximal word, provided as a string if
                `casted=False`, or as a processed (casted) word using the `cast_word`
                method from `Alphabet`.
            indices (list[int]): A list of indices that are converted into edges to construct a block.
            which (str, optional): Which facets to produce, must be one of 'upper', 'lower',
                or 'all'. Defaults to 'all'.
            casted (bool, optional): Specifies whether the word has been pre-casted using
                `Alphabet`. If False, it is treated as a string. Defaults to False.

        Returns:
            List[Block]: A list of facets.
        """

        if not casted:
            word = self._alphabet.cast_word(word)

        factors = list(word.as_coeff_mul()[1])
        upper_factors = []
        lower_factors = []

        length = len(factors)
        for i in range(length):
            base, exp = factors[i].as_base_exp()
            if i in indices:
                upper_factor = self._alphabet.get(str(base), 1)
                if exp == 1:
                    lower_factor = int(1)
                else:
                    lower_factor = base ** (exp - 1)

                upper_factors.append(lower_factor * upper_factor)
                if which in ["all", "lower"]:
                    lower_factors.append(lower_factor)
            else:
                upper_factors.append(factors[i])
                lower_factors.append(factors[i])

        upper_facets = []
        lower_facets = []
        for i in indices:
            if which in ["all", "upper"]:
                block_exp = sym.prod(
                    upper_factors[:i] + [factors[i]] + upper_factors[i + 1:],
                    start=int(1),
                )
                block = Block(str(block_exp), prod_symbol="*",
                              alphabet=self._alphabet)
                # Fix prod_symbol pylint: disable=protected-access
                block._prod_symbol = self._prod_symbol
                upper_facets.append(block)
            if which in ["all", "lower"]:
                block_exp = sym.prod(
                    upper_factors[:i] + [lower_factors[i]] +
                    upper_factors[i + 1:],
                    start=int(1),
                )
                block = Block(str(block_exp), prod_symbol="*",
                              alphabet=self._alphabet)
                lower_facets.append(block)
                # Fix prod_symbol pylint: disable=protected-access
                block._prod_symbol = self._prod_symbol

        # Filter 0's in lower facets
        if which in ["all", "lower"]:
            while Block() in lower_facets:
                lower_facets.remove(Block())

        # Return the result
        if which == "upper":
            return upper_facets
        if which == "lower":
            return lower_facets
        return upper_facets + lower_facets

    def _block_maximal_word(self, word, indices, casted=False):
        """
        Generates block based on the given maximal word and indices.

        Args:
            w (str or SymPy expression): The maximal word, provided as a string if `casted=False`,
                or as a processed (casted) word using the `cast_word` method from `Alphabet`.
            indices (list[int]): A list of indices that are converted into edges to construct a block.

        Returns:
            Block: The block represented by w(indices).
        """

        if not casted:
            word = self._alphabet.cast_word(word)

        factors = list(word.as_coeff_mul()[1])
        new_factors = []

        length = len(factors)
        for i in range(length):
            base, exp = factors[i].as_base_exp()
            if i in indices:
                upper_factor = self._alphabet.get(str(base), 1)
                if exp == 1:
                    lower_factor = int(1)
                else:
                    lower_factor = base ** (exp - 1)

                new_factors.append(lower_factor * upper_factor)
            else:
                new_factors.append(factors[i])

        block_exp = sym.prod(new_factors, start=int(1))
        block = Block(str(block_exp), prod_symbol="*", alphabet=self._alphabet)
        # Fix prod_symbol
        block._prod_symbol = self._prod_symbol  # pylint: disable=protected-access
        return block

    def _lower_facets_maximal_word_as_pairs(self, word, indices, casted=False):
        """
        Generates the lower facets of a block in word-indices representation,
        as pairs of maximal word and indices.

        Args:
            w (str or SymPy expression: The maximal word, provided as a string if `casted=False`,
                or as a processed (casted) word using the `cast_word` method from `Alphabet`.
            indices (list[int]): A list of indices that are converted into edges to construct a block.
            casted (bool, optional): Specifies whether the word has been pre-casted using `Alphabet`.
                If False, the word is treated as a string. Defaults to False.

        Returns:
            Returns:
                list[tuple[Expr, list[int]]]: List of pairs `(w, l)`, where `w` is a SymPy
                expression, and `l` is a list of indices giving the word-indices
                representation of the lower facets.
        """
        if not casted:
            word = self._alphabet.cast_word(word)

        factors = list(word.as_coeff_mul()[1])
        upper_factors = []
        lower_factors = []

        base_exp_pairs = [factor.as_base_exp() for factor in factors]
        length = len(base_exp_pairs)

        for i in range(length):
            base = base_exp_pairs[i][0]
            exp = base_exp_pairs[i][1]
            if i in indices:
                upper_factor = self._alphabet.get(str(base), 1)
                if exp == 1:
                    lower_factor = int(1)
                else:
                    lower_factor = base ** (exp - 1)

                upper_factors.append(lower_factor * upper_factor)
                lower_factors.append(lower_factor)
            else:
                upper_factors.append(factors[i])
                lower_factors.append(factors[i])
        lower_facets = []

        for i in indices:
            exp_i = base_exp_pairs[i][1]
            if i == 0:
                base_before = int(1)
            else:
                base_before = base_exp_pairs[i - 1][0]
            if i + 1 >= length:
                base_after = int(1)
            else:
                base_after = base_exp_pairs[i + 1][0]
            # Check if the resulting block is valid
            if (
                ((i - 1) in indices)
                and ((i + 1) in indices)
                and (exp_i) == 1
                and (base_before == base_after)
            ):
                continue

            # Compute maximal word
            w_max = sym.prod(
                factors[:i] + [lower_factors[i]] + factors[i + 1:], start=int(1)
            )

            # Compute the new indices
            new_indices = [
                _compute_new_index_lower_facet(
                    j, i, exp_i, base_before, base_after)
                for j in indices
                if i != j
            ]

            lower_facets.append((w_max, new_indices))

        return lower_facets

    def _check_edge(self, word1, word1_extended, word2_extended):
        """
        Internal helper function used to check whether two vertices form a 1-block
        and, if so, construct that 1-block. Assummes |word1|=|word2|+1.

        Args:
            word1 (SymPy expression or str): Base word corresponding to the first block.
                It is assumed that |w_1| = |w_2| + 1.
            word1_extended (list of tuple): List of (letter, exponent) pairs representing
                the expanded form of `word1`.
            word2_extended (list of tuple): List of (letter, exponent) pairs representing
                the expanded form of the second word.

        Returns:
            tuple:
                - bool: True if `word1` and `word2` are connected.
                - Block or None: The corresponding 1-block, if an edge is found.
                - int or None: The index `j` such that the pair (w_1, [j]) gives
                  the 1-block in word–indices form.
        """
        word1_extended_extra = [
            letter for letter, count in word1_extended for _ in range(count)
        ]
        word2_extended_extra = [
            letter for letter, count in word2_extended for _ in range(count)
        ]
        i = 0
        for j, pair in enumerate(word1_extended):
            if (
                word2_extended_extra
                == word1_extended_extra[:i] + word1_extended_extra[i + 1:]
            ):
                return True, self._block_maximal_word(word1, [j], True), j
            i += pair[1]
        return False, None, None

    def _compute_one_skeleton(
        self,
        list_words,
        list_heights,
        already_blocks,
        max_dim=1,
        alphabet=None,
        prod_symbol=None,
        verbose=False,
        check_duplicates=True,
    ):
        """
        This is an internal method used to initialize the filtration at dimensions 0 and 1.
        The filtration and incidence dictionaries are updated in place.

        Args:
            list_words (list of str or Block objects): vertices for the filtration.
            list_heights (list of int): heights associated to the words.
            already_blocks (bool): If True, assumes the words are already block objects;
                otherwise, they are converted to blocks.
            max_dim (int or None): Maximum dimension to compute. Default is 1.
            alphabet (Alphabet or None):  If None, uses `self._alphabet`.
            prod_symbol (str or None): Symbol used for concatenation in word representations.
            verbose (bool, optional): If True, prints progress messages. Default is False.
            check_duplicates (bool, optional): If True, checks that all words are unique.
                Defaults to True.
        """

        if alphabet is None:
            alphabet = self._alphabet
        elif not isinstance(alphabet, Alphabet):
            raise ValueError(
                "alphabet must be a valid Alphabet object or None")
        elif isinstance(self._alphabet, Alphabet):
            alphabet.update_letters(self._alphabet.letters_str)

        if verbose:
            print("Computing dimension 0 \u2705")
        if already_blocks:
            list_vertices = list_words
            alphabet = _combine_blocks_alphabet(list_words, self._alphabet)
        else:
            list_vertices, alphabet, prod_symbol = _convert_words_to_blocks(
                [str(w) for w in list_words],
                prod_symbol=prod_symbol,
                alphabet=self._alphabet,
            )
        self._alphabet = alphabet
        self._prod_symbol = prod_symbol

        # Check each word is provided only once
        if check_duplicates:
            if len(list_vertices) != len(set(list_vertices)):
                raise ValueError("List of words has duplicated words")
        # Initialize Dictionary at dim=0
        self.filtration_dict[0] = dict(zip(list_vertices, list_heights))

        if max_dim == 0:
            return

        if verbose:
            print("Computing dimension 1", end="")
        # Computing Graph (dim=1)
        self.filtration_dict[1] = dict()
        self._incidence_dict[1] = dict()
        list_vertices.sort(
            key=lambda x: _expand_symp_word(
                x.get_expression())[1],
            reverse=True)
        list_expanded_vertices = [
            _expand_symp_word(x.get_expression()) for x in list_vertices
        ]
        for i, word1 in enumerate(list_vertices):
            word1_extended, word1_len = list_expanded_vertices[i]
            for j in range(i + 1, len(list_vertices)):
                word2 = list_vertices[j]
                word2_extended, word2_len = list_expanded_vertices[j]
                if word1_len - word2_len == 1:
                    is_edge, block, block_index = self._check_edge(
                        word1.get_expression(), word1_extended, word2_extended
                    )
                    if is_edge:
                        self.filtration_dict[1][block] = max(
                            self.filtration_dict[0][word1],
                            self.filtration_dict[0][word2],
                        )
                        if word1.get_expression() in self._incidence_dict[1]:
                            self._incidence_dict[1][word1.get_expression()].append(
                                [block_index]
                            )
                        else:
                            self._incidence_dict[1][word1.get_expression()] = [
                                [block_index]
                            ]
                elif word1_len - word2_len > 1:
                    break
        if len(self.filtration_dict[1]) == 0:
            del self.filtration_dict[1]

    def compute_d_skeleton(
        self,
        W,
        heights=None,
        max_dim=10,
        alphabet=None,
        prod_symbol=None,
        check_duplicates=True,
        already_blocks=False,
        verbose=False,
    ):
        """
        Compute the d-skeleton of the Insertion Chain Complex generated by a set of words, C[W].
        This method replaces any existing data in the Filtration with a new complex supported on `W`.

        This method constructs all valid blocks up to the specified maximum dimension (`max_dim`)
        for a given set of words `W` with associated heights. It begins by computing the
        0- and 1-skeletons (vertices and edges), then iteratively extends to higher dimensions.

        It updates the internal `filtration_dict` to the blocks supported on `W`.

        Args:
            W (list of str or Block): List of words (or blocks, if `already_blocks=True`)
                forming the base of the complex.
            heights (list of float, optional): Height values associated with each word.
                If None, defaults to ones.
            max_dim (int, optional): Maximum dimension of the skeleton to compute. Defaults to 10.
            alphabet (Alphabet, optional): Alphabet object used together with the internal
                `self._alphabet` and any letters inferred from `W`. If provided, its symbols
                are merged with `self._alphabet`; otherwise, the new symbols are inferred entirely
                from the given words.
            prod_symbol (str, optional): Product symbol for block construction ('*', '.', or '').
                If None, inferred automatically.
            check_duplicates (bool, optional): Whether to verify that input words are unique.
                Defaults to True.
            already_blocks (bool, optional): If True, assumes the input `W` is already a list of
                `Block` objects instead of strings. Defaults to False.
            verbose (bool, optional): If True, prints progress information during computation.

        Example:
            >>> W = ['a*b', 'a*b*b', 'a*a*b','']
            >>> K = Filtration()
            >>> K.compute_d_skeleton(W, heights=[0.1, 0.3, 0.2,0.4], max_dim=2)
            >>> K[1]
            {a*b*(1,b): 0.3, a*(1,a)*b: 0.2}
        """
        if heights is None:
            heights = [1] * len(W)
        else:
            if len(heights) != len(W):
                raise ValueError(
                    "List of heighst must be same length as list of words, or None"
                )
        self.filtration_values = list(set(heights))
        self.filtration_values.sort()

        # Restart dictionaries
        self.filtration_dict = {}

        # Used when computing the d-skeleton
        self._incidence_dict = {}

        # Compute one Skeleton
        self._compute_one_skeleton(
            W,
            heights,
            already_blocks,
            max_dim,
            alphabet,
            prod_symbol,
            verbose,
            check_duplicates,
        )

        if max_dim == 1:
            self.dim = max(self.filtration_dict)
            return

        dim = 2
        while (dim <= max_dim) and len(self._incidence_dict[dim - 1]) > 0:
            if verbose:
                print(" \u2705\nComputing dimension", dim, end="")

            self._incidence_dict[dim] = dict()
            for word in self._incidence_dict[dim - 1]:
                indices = self._incidence_dict[dim - 1][word]
                possible_indices = _candidate_indices_next(indices)

                for k_indices in possible_indices:
                    lower_facets_pairs = self._lower_facets_maximal_word_as_pairs(
                        word, k_indices, True)

                    # Check if lower_facets are all there
                    all_facets_bool = True
                    for w_max, indices_face in lower_facets_pairs:
                        if (
                            w_max in self._incidence_dict[dim - 1]
                            and indices_face in self._incidence_dict[dim - 1][w_max]
                        ):
                            continue
                        else:
                            all_facets_bool = False
                            break
                    if all_facets_bool:
                        # Add to incidence dictionary
                        if word in self._incidence_dict[dim]:
                            self._incidence_dict[dim][word].append(k_indices)
                        else:
                            self._incidence_dict[dim][word] = [k_indices]

                        # Compute height of block and add
                        facets = self._facets_maximal_word(
                            word, k_indices, which="all", casted=True
                        )
                        height = max(
                            [self.filtration_dict[dim - 1][blk]
                                for blk in facets]
                        )

                        block = self._block_maximal_word(word, k_indices, True)
                        if dim in self.filtration_dict:
                            self.filtration_dict[dim][block] = height
                        else:
                            self.filtration_dict[dim] = {block: height}
            dim += 1
        self.dim = max(self.filtration_dict)

        if verbose:
            print(" \u274C")
        del self._incidence_dict

    def add_blocks(
            self,
            list_blocks,
            list_heights=None,
            prod_symbol=None,
            already_blocks=False,
            update_values=False):
        """
        Add new blocks to the Filtration.

        Extends the current Filtration by inserting additional blocks and their faces.
        This method allows dynamically modifying an existing filtration while ensuring,
        as much as possible, that the result remains a valid filtration (i.e., if
        α ≤ β, then F.get_complex(α) ⊆ F.get_complex(β)).

        Intended for expert use only, since the resulting structure may not always be
        a full Insertion Chain Complex C[W], but rather a subcomplex if some supporting
        blocks are missing. The result may depend on the order in which the blocks are
        provided.

        The behavior depends on the value of ``update_values``:

            - If ``update_values=True``, the heights of existing faces and super-faces
              are updated as needed to maintain the filtration property. This ensures
              that the new blocks can be inserted with the provided heights.

            - If ``update_values=False`` (default), the method adds the blocks only if
              their heights are consistent with the current filtration. Otherwise, the
              maximum among the heights of their faces is used instead.


        Args:
            list_blocks (list[Block] or list[str]):
                List of blocks to be added to the Filtration. If ``already_blocks`` is
                False (default), the elements are assumed to be strings representing
                blocks and will be converted. If True, they are assumed to be existing
                ``Block`` objects.

            list_heights (list[float] or float or None, optional):
                Heights assigned to each block. If a single numeric value is provided,
                it is used for all blocks. If None, all blocks receive height 1.

            prod_symbol (str or None, optional):
                Product symbol used in block representation (e.g., '*', '.', or
                ''). If not specified, it is inferred from the input blocks.

            already_blocks (bool, optional):
                If True, elements of ``list_blocks`` are assumed to be valid ``Block``
                objects. If False (default), the method attempts to convert the input
                into blocks.

            update_values (bool, optional):
                If True, existing heights of faces and super-faces are updated as needed
                to maintain consistency when inserting the new blocks. If False (default),
                only new faces are added and the provided heights are applied when valid;
                otherwise, the lowest consistent height is used instead.

        Raises:
            ValueError:
                - If ``list_heights`` is a list whose length does not match ``list_blocks``.
                - If ``list_heights`` is not a list, a numeric value, or None.
                - If the height of a block is lower than the height of one of its faces
                  already present in the filtration (when ``update_values=False``).

        Notes:
            - The internal alphabet and product symbol are updated to ensure consistency.
            - The resulting filtration may vary depending on the order in which blocks
              are inserted.
            - Updating an existing filtration in this way may be more computationally
              expensive than reconstructing a new Filtration directly from a set of words.
        """
        # Convert into blocks if needed
        if already_blocks:
            alphabet = _combine_blocks_alphabet(list_blocks, self._alphabet)
        else:
            list_blocks, alphabet, prod_symbol = _convert_words_to_blocks(
                list_blocks, prod_symbol=prod_symbol, alphabet=self._alphabet
            )
        self._alphabet = alphabet

        # Uniformalize prod_symbols pylint: disable=protected-access
        new_prods = [blk._prod_symbol for blk in list_blocks] + \
            [self._prod_symbol]
        if "*" in new_prods:
            prod_symbol = "*"
        elif "." in new_prods:
            prod_symbol = "."
        else:
            prod_symbol = ""
        self._prod_symbol = prod_symbol
        for blk in list_blocks:
            blk._prod_symbol = self._prod_symbol

        if list_heights is None:
            list_heights = [1] * len(list_blocks)
        elif isinstance(list_heights, list):
            # Check lengths agree
            if len(list_blocks) != len(list_heights):
                raise ValueError(
                    "list_heights must be same length as list_blocks")
        else:
            try:
                height = float(list_heights)
                list_heights = [height] * len(list_blocks)
            except BaseException:
                raise ValueError(
                    "list_heights must be a list, a numeric value, or None"
                )

        for block, height in zip(list_blocks, list_heights):
            dim = block.dim

            # Ensure a dictionary exists for this dimension
            if dim not in self.filtration_dict:
                self.filtration_dict[dim] = {}

            # Add or update the main block
            current_height = self.filtration_dict[dim].get(block)
            if current_height is None or update_values:
                self.filtration_dict[dim][block] = height

            # Update or add all faces
            for face in block.get_all_faces(include_self=False):
                f_dim = face.dim
                if f_dim not in self.filtration_dict:
                    self.filtration_dict[f_dim] = {}

                current_height = self.filtration_dict[f_dim].get(face)

                if current_height is None:
                    # New face — assign current block's height
                    self.filtration_dict[f_dim][face] = height
                elif update_values:
                    # Keep the lowest (earliest) height to preserve filtration
                    # order
                    self.filtration_dict[f_dim][face] = min(
                        current_height, height)

        # Recompute dimension
        self.dim = max(self.filtration_dict, default=-1)

        # Fix values so the result is a filtration
        for dim in range(1, self.dim + 1):
            for block in self.filtration_dict[dim]:
                facets = block.get_all_facets()
                all_heights = [self.filtration_dict[dim][block]] +\
                    [self.filtration_dict[dim - 1][facet] for facet in facets]
                new_height = max(all_heights)
                self.filtration_dict[dim][block] = new_height

        # Add filtration values
        self.filtration_values += list_heights
        self.filtration_values = list(set(self.filtration_values))
        self.filtration_values.sort()

    def remove_blocks(
            self,
            list_blocks,
            prod_symbol=None,
            include_upfaces=True,
            already_blocks=False):
        """
        Remove blocks from the Filtration.

        Deletes specified blocks and optionally their super-faces from the Complex.
        Intended for expert use only, since the resulting structure may not be a full
        Insertion Chain Complex C[W], but rather a subcomplex if some supported blocks
        are missing.

        Args:
            list_blocks (list of Block or string): A list of blocks to remove. If
                `already_blocks` is False (default), the elements are assumed to be strings
                representing blocks and will be converted. If True, they are assumed to be
                existing Block objects.
            prod_symbol (str or None, optional): Product symbol used in block
                representation ( '*', '.', or ''). If not specified, it is inferred
                from the input blocks.
            include_upfaces (bool, optional): If True, all super faces of the specified blocks
                are also removed, so the result is a subcomplex. Default is True.
            already_blocks (bool, optional): If True, the elements of `list_blocks`
                are assumed to be valid Block objects. If False (default), the method
                attempts to convert the input into blocks.
        """
        # Make sure it is a list
        if not isinstance(list_blocks, list):
            raise TypeError("list_blocks must be a list")

        # Convert into blocks if needed
        if not already_blocks:
            list_blocks, _alphabet, prod_symbol = _convert_words_to_blocks(
                list_blocks, prod_symbol=prod_symbol, alphabet=self._alphabet
            )

        # Dictionary of blocks to remove
        blocks_to_remove = {i: [] for i in range(self.dim + 1)}
        for block in list_blocks:
            if block.dim in self.filtration_dict and block in self.filtration_dict[block.dim]:
                blocks_to_remove[block.dim].append(block)

        # Find super-faces if needed
        if include_upfaces:
            for dimension in range(1, self.dim + 1):
                for block in self.filtration_dict[dimension]:
                    facets = block.get_all_facets()
                    if any(
                            facet in blocks_to_remove[dimension - 1] for facet in facets):
                        blocks_to_remove[dimension].append(block)

        # Remove blocks
        for dimension in blocks_to_remove:
            for block in blocks_to_remove[dimension]:
                if block in self.filtration_dict[dimension]:
                    del self.filtration_dict[dimension][block]

        # Update filtration_dict
        for dimension in range(self.dim, -1, -1):
            if len(self.filtration_dict[dimension]) == 0:
                del self.filtration_dict[dimension]
            else:
                break

        # Update filtration values
        self.filtration_values = list({
            val
            for dimension_dict in self.filtration_dict.values()
            for val in dimension_dict.values()
        })
        self.filtration_values.sort()

        # Recompute dimension
        self.dim = max(self.filtration_dict, default=-1)

    def get_complex(self, height=None, max_dim=None):
        """
        Constructs and returns a Complex object that includes all blocks
        from the filtration whose height is less than or equal to the specified
        value. The construction can also be limited to a specified maximum dimension.

        Args:
            height (float or int, optional): The maximum filtration value to include.
                If None, the largest available filtration value is used.
            max_dim (int, optional): The maximum dimension to include in the complex.
                If None, the full dimension of the filtration is used.

        Returns:
            Complex: A Complex object containing all blocks up to the specified
                height and dimension.
        """
        from InDelsTopo.complex import Complex  # pylint: disable=import-outside-toplevel

        if height is None:
            height = max(self.filtration_values, default=np.inf)
        if max_dim is None:
            max_dim = self.dim
        max_dim = min(self.dim, max_dim)

        complex_dict = {
            dim: [
                block
                for block in self.filtration_dict[dim]
                if self.filtration_dict[dim][block] <= height
            ]
            for dim in range(max_dim + 1)
        }
        for dim in range(max_dim, -1, -1):
            if len(complex_dict[dim]) == 0:
                del complex_dict[dim]
            else:
                break
        return Complex(
            alphabet=self._alphabet,
            prod_symbol=self._prod_symbol,
            complex_dict=complex_dict,
            height=height,
        )

    def get_euler_curve(self, x_values=None):
        """
        This method evaluates the Euler characteristic of the complex at different
        filtration heights and returns the resulting curve as paired x- and y-values.

        Args:
            x_values (list of float or int, optional): Filtration heights at which to
                compute the Euler characteristic. If None, all existing filtration
                values are used.

        Returns:
            tuple:
                - list of float or int: Sorted filtration heights (x-values).
                - list of int: Corresponding Euler characteristic values (y-values).
        """
        if x_values is None:
            x_values = self.filtration_values
        x_values.sort()
        y_values = [self.get_complex(h).euler_characteristic()
                    for h in x_values]
        return x_values, y_values

    def get_persistent_homology_barcodes(
        self, max_dim=None, inf_value=np.inf, get_height_indices=False
    ):
        """
        Compute persistent homology barcodes for the filtration using Z2 coefficients.

        The method performs a lower-star filtration on the complex according to the
        vertex heights and computes persistent homology up to the specified dimension.

        Args:
            max_dim (int, optional): Maximum dimension to compute accurately.
                The skeleton of dimension up to max_dim+1 is used for this computation
                if it was previously computed. If None or greater than the filtration's
                dimension, all dimensions are included.
            inf_value (float, optional): Value to assign to features that do not die
                within the filtration. Defaults to infinity.
            get_height_indices (bool, optional): If True, also return the indices
                corresponding to the birth and death heights in the filtration. Defaults to False.

        Returns:
            dict or tuple:
                - If `get_height_indices=False`: A dictionary mapping dimension `d`
                  to a list of tuples `(birth, death)` representing persistent homology intervals.
                - If `get_height_indices=True`: A tuple of two dictionaries:
                    1. Barcodes as above.
                    2. Corresponding indices of birth and death heights in `filtration_values`.

        Notes:
            The algorithm uses a boundary matrix over Z2 and reduces it following
            the standard persistence algorithm (see arxiv:1506.08903).
            Features that never die are assigned `inf_value` as their death time.
        """
        # Import csc_matrix
        from scipy.sparse import csc_matrix  # pylint: disable=import-outside-toplevel

        if max_dim is None or max_dim >= self.dim:
            max_dim = self.dim
            using_skeleton = False
        elif max_dim == self.dim - 1:
            max_dim = self.dim
            using_skeleton = True
        else:
            max_dim += 1
            using_skeleton = True

        # Order the blocks according to their dimension and height
        ordered_blocks = []
        for dimension in range(max_dim + 1):
            ordered_blocks += self[dimension]
        ordered_blocks.sort(key=lambda blk: (self[blk.dim][blk], blk.dim))

        ordered_blocks_dict = {
            ordered_blocks[i]: i for i in range(len(ordered_blocks))}
        heights = [self[blk.dim][blk] for blk in ordered_blocks]

        # Construct the boundary matrix
        cols = []
        rows = []
        data = []
        for id_col, blk in enumerate(ordered_blocks):
            facets = blk.get_all_facets()
            for facet in facets:
                cols.append(id_col)
                rows.append(ordered_blocks_dict[facet])
                data.append(1)

        N = len(ordered_blocks)
        boundary_matrix = csc_matrix(
            (data, (rows, cols)), shape=(N, N)).tolil()

        # Perform row reduction (we follow the algorithm from
        # https://arxiv.org/pdf/1506.08903)
        low = []  # maps row index to pivot column index
        barcodes = {}

        if get_height_indices:
            barcodes_indices = {}
            height_indices = {self.filtration_values[i]: i for i in range(
                len(self.filtration_values))}

        for j in range(N):
            col = boundary_matrix.getcol(j)
            if col.nnz > 0:
                # Set initial value for low(j)
                low_j = col.nonzero()[0][-1]
                while low_j in low:
                    i = low.index(low_j)
                    # add column i to column j
                    col = col + boundary_matrix.getcol(i)
                    col.data = col.data % 2  # Make adition modulo 2
                    col.eliminate_zeros()
                    if col.nnz > 0:
                        low_j = col.nonzero()[0][-1]  # update low_j value
                    else:
                        low_j = -1
                        break
                boundary_matrix[:, j] = col  # update column j in the matrix
                low.append(low_j)  # Save value for low_j
            else:
                # Set -1, for undefined low(j)
                low.append(-1)

        # Extract (birth, death) pairs
        for j, low_j in enumerate(low):
            if low_j >= 0:
                birth = heights[low_j]
                death = heights[j]
                if death > birth:
                    dim = ordered_blocks[low_j].dim
                    if dim not in barcodes:
                        barcodes[dim] = []
                        if get_height_indices:
                            barcodes_indices[dim] = []
                    barcodes[dim].append((birth, death))
                    if get_height_indices:
                        barcodes_indices[dim].append(
                            (height_indices[birth], height_indices[death])
                        )
            elif j not in low:
                birth = heights[j]
                dim = ordered_blocks[j].dim
                if dim not in barcodes:
                    barcodes[dim] = []
                    if get_height_indices:
                        barcodes_indices[dim] = []
                barcodes[dim].append((birth, inf_value))
                if get_height_indices:
                    barcodes_indices[dim].append(
                        (height_indices[birth], np.inf))

        if using_skeleton:
            if max_dim in barcodes:
                del barcodes[max_dim]
                if get_height_indices:
                    del barcodes_indices[max_dim]

        if get_height_indices:
            return barcodes, barcodes_indices
        return barcodes

    def get_graph(
        self,
        height=None,
        height_id=None,
        show_labels=True,
        max_dim=5,
        positions=None,
        initial_positions=None,
        fixed=None,
        recompute=False,
        colors_by_dim=None,
        ax=None
    ):
        """
        Generate a graphical representation of the filtration up to a specified dimension.

        Positions of vertices can be computed automatically or provided manually.
        Only accurate for low-dimensional complexes (typically dim <= 3).

        Args:
            height (float, optional): Maximum height value for including blocks in the graph.
                Defaults to the maximum filtration value.
            height_id (int, optional): Index into the sorted list of filtration values.
                Used if `height` is None. Defaults to None.
            show_labels (bool, optional): Whether to display labels on the vertices. Defaults to True.
            max_dim (int, optional): Maximum dimension of blocks to include in the graph. Defaults to 5.
            positions (dict, optional): Dictionary of vertex positions.
                If None, positions are computed automatically. Once computed, they are reused
                everytime this method is called, unless recompute is set to True.
            initial_positions (dict, optional): Initial positions used to seed the
                automatic layout algorithm.
            fixed (list or None, optional): List of vertex keys to fix in place when computing positions.
                Defaults to None.
            recompute (bool, optional): Whether to recompute vertex positions even
                if already stored. Defaults to False.
            colors_by_dim (list of str, optional): List of colors to use for each dimension.
                If None, defaults to ['black', 'gray', 'yellow', 'red', 'blue', 'purple'].
            ax (matplotlib.axes._subplots.Axes3DSubplot, optional): A Matplotlib Axes
                object to draw the plot on. If None, a new figure and axes are created.
                Defaults to None.

        Returns:
            matplotlib.axes.Axes: Matplotlib axes object containing the drawn graph.
        """
        if self.dim == -1:
            return None
        if (positions is None) or recompute:
            if (self._positions_dict is None) or recompute:
                self._positions_dict = graphics.compute_vertex_positions(
                    self, pos0=initial_positions, fixed=fixed
                )
            # Make sure position includes all vertices
            elif any(vertex not in self._positions_dict for vertex in self[0]):
                # Update with pos0
                if isinstance(initial_positions, dict):
                    for vertex in initial_positions:
                        self._positions_dict[vertex] = initial_positions[vertex]
                # Compute for all vertices
                self._positions_dict = graphics.compute_vertex_positions(
                    self, pos0=self._positions_dict, fixed=list(
                        self._positions_dict.keys())
                )

            # Get the positions
            positions = self._positions_dict

        if height is None:
            try:
                height = self.filtration_values[height_id]
            except BaseException:
                height = max(self.filtration_values, default=np.inf)

        ax = graphics.make_graph(
            self,
            pos=positions,
            show_labels=show_labels,
            max_dim=max_dim,
            height=height,
            already_complex=False,
            colors_by_dim=colors_by_dim,
            ax=ax
        )

        return ax

    def __str__(self):
        to_print = "Filtration of Insertion Chain Complexes: \n"
        to_print += "alphabet: " + str(self._alphabet) + ".\n"
        try:
            filtration_values_range = (
                "["
                + str(self.filtration_values[0])
                + ","
                + str(self.filtration_values[-1])
                + "]"
            )
        except BaseException:
            filtration_values_range = "[]"
        to_print += "heights in: " + filtration_values_range + ".\n"
        to_print += "dimension: " + str(self.dim) + ".\n"
        to_print += "vertices: " + str(len(self[0])) + ".\n"
        to_print += (
            "blocks: " + str(sum([len(self[k])
                             for k in range(self.dim + 1)])) + "."
        )
        return to_print

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        if 0 <= key <= self.dim:
            return self.filtration_dict[key]
        return {}

    def get_alphabet(self):
        """Returns the alphabet attribute. """
        return self._alphabet

    def get_block_height(
            self,
            block,
            already_block=False):
        """
        Returns the height of the block in the filtration or None if it is not in the filtration. 

        Args:
            block (Block or string): A block to check its height. If `already_blocks` is False (default), 
            it is assumed to be a string representing the block. If True, it is assumed to be an
                existing Block object.
            already_block (bool, optional): If True, `block`
                is assumed to be a valid Block object. If False (default), the method
                attempts to convert the input into a block.
        """
        # Convert into blocks if needed
        if not already_block:
            block=Block(block, prod_symbol=self._prod_symbol, alphabet=self._alphabet)
        
        dim=block.dim
        if block in self[dim]:
            return self[dim][block]
        return None
