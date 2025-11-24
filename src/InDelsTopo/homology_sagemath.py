"""
Module to compute homology with integer coefficients for Insertion Chain Complexes.

This module requires SageMath to run and cannot be executed in a standard
Python environment.
"""

from sage.rings.integer_ring import ZZ                # pylint: disable=import-error
from sage.homology.chain_complex import ChainComplex  # pylint: disable=import-error
from sage.matrix.constructor import matrix            # pylint: disable=import-error

import numpy as np


def sort_expressions(blocks_list):
    """
    Sort a list of blocks by the string representation of their expressions.

    Args:
        blocks_list (list[Block]): List of Blocks to sort.

    Returns:
        list[Block]: The sorted list of Blocks.
    """
    blocks_list.sort(key=lambda x: str(x.get_expression()))
    return blocks_list


def create_chain_complex(blocks, get_ordered_blocks=False):
    """
    Create a SageMath chain complex from a dictionary of blocks.

    This function builds boundary matrices for each dimension of the insertion
    chain complex and returns a ChainComplex object. Optionally, it can also
    return the blocks sorted within each dimension.

    Args:
        blocks (dict[int, iterable[Block]]): A dictionary mapping dimension `d`
            to an iterable of Block objects of that dimension.
        get_ordered_blocks (bool, optional): If True, also return the blocks
            ordered by their expressions within each dimension. Defaults to False.

    Returns:
        ChainComplex or tuple[ChainComplex, dict[int, list[Block]]]:
            - If `get_ordered_blocks` is False, returns a ChainComplex object
              representing the boundary operators.
            - If `get_ordered_blocks` is True, returns a tuple `(ChainComplex, blocks_ordered)`,
              where `blocks_ordered` is a dictionary of lists of blocks sorted by expression.
    """
    blocks_ordered = {d: sort_expressions(list(blocks[d])) for d in blocks}
    boundary_operators = [matrix(ZZ, np.ones((1, len(blocks_ordered[0]))))]
    for d in range(1, max(blocks) + 1):
        matrix_boundary = matrix(
            ZZ, np.zeros((len(blocks_ordered[d - 1]), len(blocks_ordered[d])))
        )
        all_faces = blocks_ordered[d - 1]
        for j, block in enumerate(blocks_ordered[d]):
            block_boundary = block.boundary()
            for term in block_boundary.get_blocks():
                i = all_faces.index(term)
                boundary_dict_blocks = block_boundary.get_dict_blocks()
                matrix_boundary[i, j] = int(boundary_dict_blocks[term])
        boundary_operators.append(matrix_boundary)
    if get_ordered_blocks:
        return (
            ChainComplex(dict(enumerate(boundary_operators)), degree=-1),
            blocks_ordered,
        )
    return ChainComplex(dict(enumerate(boundary_operators)), degree=-1)
