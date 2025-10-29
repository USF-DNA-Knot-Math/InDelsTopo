# InDelsTopo 
A Python package to analyze topological properties of sets of words when their main source of variation are insertions and deletions, using the Insertion Chain Complex framework.

## Introduction

**InDelsTopo** implements the **Insertion Chain Complex** framework (Jonoska, Martinez-Figueroa, and Saito, [arXiv:2509.12607](https://arxiv.org/abs/2509.12607)), originally developed to analyze variation in collections of DNA sequences following double-strand break repair. 

The package enables general analysis of word sets where **insertions and deletions** are the primary source of variation, using topological tools such as **persistent homology** and **Euler characteristic curves**. 

---

## Features

- Define and generate Insertion Chain Complexes and Filtrations based on sets of words, over any finite alphabet. 
- Support for computing homology using either:
  - Pure Python/NumPy backend (with coefficients in $\mathbb{Z}_2$). 
  - [SageMath](https://www.sagemath.org/) for symbolic computation with coefficients in any supported ring (including $\mathbb{Z}$).
- Computation of persistent homology barcodes for Filtrations.
- Computation of Euler Characteristic Curves for Filtrations.
- Graphic representation of Complexes (mostly useful for small and simple complexes). 
- Functions and classes to create and manipulate Blocks.

---

## Requirements

- Python 3.8+
- NumPy
- Networkx
- Simpy
- Scipy

**Optional (for extended functionality):**
- [SageMath](https://www.sagemath.org/) (for symbolic homology computation)


## Installation

```bash
pip install InDelsTopo
# Required dependencies (if not already installed):
# pip install sympy numpy matplotlib networkx scipy
```

## Quick Example 

```python
from InDelsTopo import Complex

# Define a set of words
W = ['abbc', 'ab^2c^2', 'abc', 'abc**2', 'abc**2d']

# Create a complex
K = Complex()
K.compute_d_skeleton(W)

# Inspect
print(K)
print("2-blocks:", K[2])

# Visualize (uncomment for interactive plot)
# %matplotlib widget
K.get_graph()
```

## Cheat Sheet: Main Classes and Methods

| Class / Object | Method / Attribute | Description |
|----------------|------------------|-------------|
| **Complex** | `compute_d_skeleton(W, max_dim=d, verbose=False)` | Construct the d-skeleton of the insertion chain complex from a set of words `W`. |
|  | `add_blocks(list_blocks, already_blocks=False)` | Add new blocks to the complex. If `already_blocks=True`, the elements are `Block` objects. |
|  | `remove_blocks(list_blocks)` | Remove specified blocks and all super-faces containing them. |
|  | `get_graph()` | Visualize the complex. |
|  | `dim` | Maximum dimension of the complex. |
|  | `euler_characteristic` | Euler characteristic of the complex. |
|  | `get_betti_numbers_z2()` | Compute Betti numbers with coefficients in $\mathbb{Z}_2$. |
|  | `get_maximal_blocks()` | Returns dictionary of maximal blocks by dimension. |
| **Filtration** | `compute_d_skeleton(W, heights=None, max_dim=d)` | Construct a filtration using word heights. |
|  | `add_blocks(list_blocks, list_heights)` | Add blocks to the filtration, optionally updating heights of faces. |
|  | `remove_blocks(list_blocks)` | Remove blocks and all higher-dimensional blocks containing them. |
|  | `get_complex(height=None, max_dim=d)` | Return a subcomplex at a given height. |
|  | `filtration_dict` | Dictionary of blocks by dimension and their heights. |
|  | `dim` | Maximum dimension of the filtration. |
|  | `filtration_values` | Sorted list of all heights in the filtration. |
| **Block** | Constructor: `Block(string)` | Create a block from a string representation. |
|  | `get_all_faces()` | Returns all subfaces of the block. |
|  | `get_all_facets()` | Returns codimension-1 faces. |
|  | `get_upper_facets()` / `get_lower_facets()` | Returns only upper or lower facets. |
|  | `get_all_vertices()` | Returns all vertices. |
|  | `get_face(indices_plus, indices_minus)` | Returns a specific face by upper/lower indices. |
|  | `get_vertex(indices)` | Returns a specific vertex by index list. |
| **Chain** | Constructor: `Chain(string)` | Linear combination of blocks with integer coefficients. |
|  | Arithmetic: `+`, `-`, `*` | Combine chains and multiply by integers. |
|  | `get_blocks()` | Returns list of blocks in the chain. |
|  | `get_dict_blocks()` | Returns dictionary of blocks with their coefficients. |
|  | `boundary()` | Compute boundary of the chain. |

## Documentation

Full documenation can be found at [https://usf-dna-knot-math.github.io/InDelsTopo/]

## Tutorial
The full Jupyter Notebook tutorial is available in [tutorials/InDelsTopo_Tutorial.ipynb]
It covers:
1. Complexes: constructing, inspecting, and visualizing.
2. Filtrations: heights, sublevel sets, and persistent homology.
3. Blocks and Chains: creating, combining, and computing boundaries.
4. Manipulating complexes and filtrations: adding/removing blocks.
5. Optional SageMath integration for homology over $\mathbb{Z}$.

## Acknowledgements 
This package was built under auspices of the Southeast Center for Mathematics and Biology, an NSF-Simons Research Center for Mathematics of Complex Biological Systems, under National Science Foundation Grant No. DMS-1764406 and Simons Foundation Grant No. 594594 as well as NSF DMS-2054321, CCF-2107267, CCF-2505771 and the W.M. Keck Foundation.  

## Citation
If you use **InDelsTopo** in your research, please cite the following preprint:

Jonoska, N., Martinez-Figueroa, F., & Saito, M. (2025). *The Insertion Chain Complex: A Topological Approach to the Structure of Word Sets*. [arXiv:2509.12607](https://arxiv.org/abs/2509.12607)

## Contributing
We welcome contributions to **InDelsTopo**! You can help by:

- Reporting bugs or issues via the [GitHub Issues](https://github.com/USF-DNA-Knot-Math/InDelsTopo/issues) page.
- Submitting pull requests with improvements, bug fixes, or new examples.
- Suggesting new features or documentation improvements.

Before contributing, please ensure that your code follows the existing style and passes tests.

## License
**InDelsTopo** is released under the MIT License. See the [LICENSE](LICENSE.txt) file for details.
