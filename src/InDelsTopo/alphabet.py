"""Module for the Alphabet class and related helper functions."""

# Third party imports
import sympy as sym


# Alphabet Class
class Alphabet:
    """
    Represents an alphabet of symbols used to build words.

    Each symbol in the alphabet is represented as a noncommutative SymPy symbol.
    For every letter, an associated "edge" symbol (or 1-block symbol) is also
    created. For example, the letter 'a' has a corresponding edge symbol '(1,a)'.
    Internally, these symbols are handled as elements of a noncommutative SymPy
    monoid, allowing symbolic manipulations of words as products.

    This class provides methods to:
    - Add individual letters.
    - Update the alphabet with multiple new letters.
    - Convert words (sequences of letters) into symbolic products of letters.

    Attributes:
        letters (dict): Maps letter strings to their SymPy symbol representation.
        edges (dict): Maps letter strings to their corresponding edge (1-block) symbols.
        letters_str (list): List of letter strings currently in the alphabet.
    """

    def __init__(self, letters_str=None):
        """
        Initialize an Alphabet with optional letters.

        Args:
            letters_str (list of str, optional): Initial letters for the alphabet.
                If None, an empty alphabet is created. Duplicate letters are removed.

        Side Effects:
            Initializes the following attributes:
                - letters (dict): Maps letter strings to SymPy symbols.
                - edges (dict): Maps letter strings to corresponding edge symbols.
                - letters_str (list): List of unique letter strings in the alphabet.
        """
        letters = {}
        edges = {}
        letters_str = list(set(letters_str)) if letters_str is not None else []
        for symbol in letters_str:
            letters[symbol] = sym.Symbol(symbol, commutative=False)
            edges[symbol] = sym.Symbol("(1,_)".replace("_", symbol), commutative=False)

        self.letters = letters
        self._edges = edges
        self.letters_str = letters_str

    def get(self, letter, dim=0):
        """
        Return the symbolic representation of a letter or its edge.

        Args:
            letter (str): The letter to retrieve.
            dim (int, optional): Dimension; 0 for the letter symbol,
                1 for the edge symbol. Defaults to 0.

        Returns:
            SymPy.Symbol: The corresponding SymPy symbol.
        """
        if dim == 0:
            return self.letters[letter]
        if dim == 1:
            return self._edges[letter]
        raise ValueError("Dimension dim must be 0 or 1.")

    def add_letter(self, symbol):
        """
        Add a new letter to the alphabet if it does not already exist.

        Args:
            symbol (str): The new letter to add.
        """
        if symbol not in self.letters_str:
            self.letters_str.append(symbol)
            self.letters[symbol] = sym.Symbol(symbol, commutative=False)
            self._edges[symbol] = sym.Symbol(
                "(1,_)".replace("_", symbol), commutative=False
            )

    def update_letters(self, letters_str):
        """
        Add multiple letters to the alphabet, ignoring duplicates.

        Args:
            letters_str (iterable of str): letters to add.
        """
        new_letters = set(letters_str).difference(self.letters_str)
        for symbol in new_letters:
            self.letters_str.append(symbol)
            self.letters[symbol] = sym.Symbol(symbol, commutative=False)
            self._edges[symbol] = sym.Symbol(
                "(1,_)".replace("_", symbol), commutative=False
            )

    def cast_word(self, word, check_letters=True):
        """
        Convert a word (sequence of letters) into a symbolic product.

        Args:
            word (iterable of str): The word to convert, as a list of symbols or a string.
            check_letters (bool, optional): If True, updates the alphabet with
                any new letters found in the word. Defaults to True.

        Returns:
            sympy.Expr: SymPy expression of the word, using the symbols in this alphabet.
        """
        if check_letters:
            self.update_letters(word)
        return sym.prod([self.letters[symbol] for symbol in word], start=int(1))

    def __str__(self):
        return "Alphabet with letters: " + str(sorted(self.letters_str))

    def __repr__(self):
        return self.__str__()
