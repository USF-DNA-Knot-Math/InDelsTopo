"""
Utility functions for the InDelsTopo package.

This module provides auxiliary functions used various classes,
including operations on alphabets, blocks, and chains that do not naturally belong
to a single class.
"""

# Standard library imports
import itertools

# Third party imports
import sympy as sym

# Local imports
from InDelsTopo.alphabet import Alphabet


# Functions
# ----------------------------------------------------------------------------------------
def _expand_symp_word(word):
    """
    Expand a SymPy word into its base-exponent pairs and total length.

    Args:
        word (str or SymPy expression): The word or expression to expand.

    Returns:
        tuple:
            list of tuple: Each element is a (base, exponent) pair.
            int: The total length, i.e., the sum of all exponents.
    """
    pairs = [pair.as_base_exp() for pair in sym.sympify(word).as_coeff_mul()[1]]
    length = sum([pair[1] for pair in pairs])
    return pairs, length


def _combine_blocks_alphabet(blocks_list, alphabet=None):
    """
    Combine the alphabets of a collection of blocks into a single Alphabet.

    This function aggregates all unique symbols in the alphabets of the blocks in blocks_list,
    merges them with an optional initial alphabet, and returns a new
    `Alphabet` object containing the combined symbols.

    Args:
        blocks_list (iterable): A collection of Block objects.
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


def _powerset(iterable):
    """
    Return the powerset of the given iterable. Uses the package itertools.
    """
    iterable = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(iterable, r) for r in range(len(iterable) + 1)
    )
