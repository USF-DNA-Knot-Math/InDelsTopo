# InDelsTopo/__init__.py
"""
InDelsTopo: A Python package to analyze topological properties of sets of words when
their main source of variation are insertions and deletions,
using the Insertion Chain Complex framework.
"""

# Version
__version__ = "0.1.0"

# Expose main classes/functions
from .block import Block
from .chain import Chain
from .alphabet import Alphabet
from .filtration import Filtration
from .complex import Complex
