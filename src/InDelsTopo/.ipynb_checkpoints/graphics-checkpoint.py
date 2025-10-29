"""Module to produce graphic representation of Insertion Chain Complexes"""

# pylint: disable=invalid-name
# Import external packages
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Auxiliary functions
def fibonacci_sphere(n):
    """Auxiliary function to produce `n` points that look `equally spaced` on a unit sphere"""
    points = []
    phi = np.pi * (3.0 - np.sqrt(5))
    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2
        radius = np.sqrt(1 - y * y)
        theta = phi * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))
    return np.array(points)


def platonic_solid(n):
    """Auxiliary function to produce the vertices of a platonic solid
    with `n` vertices inscribed on a unit sphere"""
    if n == 4:  # Tetrahedron
        return np.array(
            [
                [1, 1, 1],
                [-1, -1, 1],
                [-1, 1, -1],
                [1, -1, -1],
            ]
        ) / np.sqrt(3)
    if n == 6:  # Octahedron
        return np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        )
    if n == 8:  # Cube
        return np.array(
            [
                [1, 1, 1],
                [1, 1, -1],
                [1, -1, 1],
                [1, -1, -1],
                [-1, 1, 1],
                [-1, 1, -1],
                [-1, -1, 1],
                [-1, -1, -1],
            ]
        ) / np.sqrt(3)
    if n == 12:  # Icosahedron
        phi = (1 + np.sqrt(5)) / 2
        verts = []
        for i in [-1, 1]:
            for j in [-1, 1]:
                verts.append([0, i * phi, j])
                verts.append([i, 0, j * phi])
                verts.append([i * phi, j, 0])
        return np.array(verts) / np.linalg.norm([phi, 1, 0])
    if n == 20:  # Dodecahedron
        phi = (1 + np.sqrt(5)) / 2
        verts = []
        for a in [-1, 1]:
            for b in [-1, 1]:
                verts.append([0, a / phi, b * phi])
                verts.append([a / phi, b * phi, 0])
                verts.append([a * phi, 0, b / phi])
        for i in [-1, 1]:
            verts.append([i, i, i])
            verts.append([i, i, -i])
            verts.append([i, -i, i])
            verts.append([-i, i, i])
        return np.array(verts) / np.linalg.norm([phi, 1, 0])
    return None


def get_points_on_sphere(n):
    """Auxiliary function to produce `n` points that look `equally spaced`
    on the unit sphere. Uses the vertices of the closest platonic solid,
    or if `n`>20, it uses the fibonacci_sphere function.
    """
    platonic_ns = [4, 6, 8, 12, 20]
    if n <= 20:
        for m in platonic_ns:
            if n <= m:
                break
        platonic_vertices = platonic_solid(m)
        return platonic_vertices[:n]
    return fibonacci_sphere(n)


def format_vertex_text(vertex):
    """Formats the vertex with LaTeX for showing it in the graph."""
    if str(vertex) == "1":
        return "$1$"
    text = ""
    factors = vertex.as_coeff_mul()[1]
    for factor in factors:
        base_pair = factor.as_base_exp()
        text += str(base_pair[0])
        if base_pair[1] > 1:
            text += "^{" + str(base_pair[1]) + "}"
    return "$" + text + "$"


def recenter(vertices, pos, prev_bound):
    """
    Recenter a set of points around a new centroid.

    Given a list of keys `vertices` and a dictionary `pos` mapping keys to coordinates,
    this function computes the centroid of the points corresponding to `vertices`,
    shifts all those points so that the centroid is moved relative to `prev_bound`
    plus a radius offset, and returns the updated positions along with the
    position of the farthest point from the centroid.

    Args:
        vertices (Iterable): Keys identifying the subset of points to recenter.
        pos (dict): Dictionary mapping keys to NumPy arrays or coordinate lists.
        prev_bound (numpy.ndarray): Reference position to shift the centroid with respect to.

    Returns:
        tuple:
            dict: Updated `pos` dictionary with recentered coordinates for keys in `vertices`.
            numpy.ndarray: Coordinates of the point in `vertices`
                that was farthest from the centroid.
    """
    vertices = list(vertices)
    coords = [pos[k] for k in vertices]
    centroid = np.mean(coords, axis=0)
    r = np.sqrt(np.max(np.sum((np.array(coords) - centroid) ** 2, axis=1)))
    i = np.argmax(np.sum((np.array(coords) - centroid) ** 2, axis=1))
    pos.update({k: pos[k] - centroid + prev_bound + r + 1 for k in vertices})
    return pos, pos[vertices[i]]


def get_starting_position(letters, vertices):
    """
    Compute initial 3D positions for a set of vertices based on letter positions.

    Each letter is assigned a point on a unit sphere, and each vertex in `vertices`
    is represented as a linear combination of letter positions according to its
    coefficients. Small random noise is added for vertices with multiple factors
    to avoid exact overlaps. The special vertex `1` (if present) is positioned at
    the origin.

    Args:
        letters (list of str): List of letters in the alphabet.
        vertices (iterable): Iterable of SymPy expressions representing vertices.

    Returns:
        dict:
            Mapping from each vertex in `vertices` to a numpy array of shape (3,)
            representing its 3D coordinates.
    """
    positions_letters = get_points_on_sphere(len(letters))
    pos_letters = {letters[i]: positions_letters[i]
                   for i in range(len(letters))}
    pos0 = {}
    for vertex in vertices:
        coord = 0
        factors = vertex.as_coeff_mul()[1]
        for factor in factors:
            base_pair = factor.as_base_exp()
            coord += pos_letters[base_pair[0]] * int(base_pair[1])
        if len(factors) > 1:
            coord = coord + np.random.normal(0, scale=0.1, size=(3,))
        pos0[vertex] = coord
    # Position for '1'
    if sym.sympify(1) in vertices:
        pos0[sym.sympify(1)] = np.array([0, 0, 0])
    return pos0


# Make the graph
def make_graph(
    K,
    pos,
    show_labels=False,
    max_dim=5,
    height=1,
    already_complex=False,
    colors_by_dim=None,
    ax=None
):
    """
    Create a 3D plot of an insertion chain complex using matplotlib.

    This function plots vertices, edges, and 2D faces of the complex in 3D space
    according to the positions provided in `pos`. Colors are assigned by dimension
    for visualization, and optional vertex labels can be displayed.

    Args:
        K (Complex or Filtration): The object representing the complex to plot.
            If `already_complex` is True, `K` is assumed to be a Complex.
            Otherwise, `K` is a Filtration from which a lower-star subcomplex is extracted.
        pos (dict): Mapping from vertices (SymPy expressions) to 3D coordinates
            as numpy arrays of shape (3,).
        show_labels (bool, optional): Whether to display vertex labels in the plot.
            Defaults to False.
        max_dim (int, optional): Maximum dimension of faces to include in the plot.
            Defaults to 5.
        height (float, optional): Height parameter used to obtain the lower-star
            subcomplex when `already_complex` is False. Defaults to 1.
        already_complex (bool, optional): If True, `K` is assumed to already be
            a Complex; otherwise, a lower-star subcomplex is generated via `get_complex`.
            Defaults to False.
        colors_by_dim (list of str, optional): List of colors to use for each dimension.
            If None, defaults to ['black', 'gray', 'yellow', 'red', 'blue', 'purple'].
        ax (matplotlib.axes._subplots.Axes3DSubplot, optional): A Matplotlib Axes
            object to draw the plot on. If None, a new figure and axes are created.
            Defaults to None. 

    Returns:
        matplotlib.axes._subplots.Axes3DSubplot:
            The 3D axes object containing the plotted complex. It can be further
            customized or displayed using `plt.show()`.
    """
    if colors_by_dim is None:
        colors_by_dim = ["black", "gray", "yellow", "red", "blue", "purple"]

    # Now plotting
    if already_complex:
        sub_complex = K.get_complex(max_dim=max_dim)
    else:
        sub_complex = K.get_complex(height=height, max_dim=max_dim)
    maximal_blocks = sub_complex.get_maximal_blocks()

    if sub_complex.dim >= 0:
        vertices = [v.get_expression() for v in sub_complex[0]]
    else:
        vertices = []

    if sub_complex.dim >= 1:
        edges = sub_complex[1]
        edges = [e.get_all_vertices() for e in edges]
    else:
        edges = []
        
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

    pos = {v: pos[v] for v in vertices}

    # Draw vertices
    for vertex in vertices:
        (x, y, z) = pos[vertex]
        ax.scatter(x, y, z, s=10, c=colors_by_dim[0])
        if show_labels:
            ax.text(
                x, y, z, format_vertex_text(vertex), fontsize=10, va="bottom", ha="left"
            )

    # Draw edges
    for u, v in edges:
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        z = [pos[u][2], pos[v][2]]
        ax.plot(x, y, z, color=colors_by_dim[1])

    # Squares
    for dim in maximal_blocks:
        square_parts = []
        if dim >= 2:
            for block in maximal_blocks[dim]:
                new_squares = [
                    s for s in block.get_all_faces(True) if s.dim == 2]
                new_squares_coords = [
                    [
                        pos[s.get_vertex([])],
                        pos[s.get_vertex([1])],
                        pos[s.get_vertex([1, 2])],
                        pos[s.get_vertex([2])],
                    ]
                    for s in new_squares
                ]
                for s in new_squares_coords:
                    centroid = np.mean(s, axis=0)
                    square_parts += [
                        [s[0], s[1], centroid],
                        [s[1], s[2], centroid],
                        [s[2], s[3], centroid],
                        [s[3], s[0], centroid],
                    ]

            # Draw squares
            poly3d = Poly3DCollection(
                square_parts, facecolors=colors_by_dim[dim], edgecolors=None, alpha=0.2
            )
            ax.add_collection3d(poly3d)

    # Set plot limits and aesthetics
    ax.set_box_aspect([1, 1, 1])  # equal aspect ratio
    ax.axis("off")  # remove axes

    return ax


# Compute initial positions of vertices
def compute_vertex_positions(K, pos0=None, fixed=None):
    """
    Compute 3D coordinates for the vertices of an insertion chain complex.

    This function assigns 3D positions to the vertices of a given complex to
    enable graphical visualization. Some vertex positions can be preset and
    kept fixed throughout the computation.

    Args:
        K (Complex):
            The insertion chain complex whose vertex coordinates are to be computed.
        pos0 (dict):
            A mapping from vertices (as SymPy expressions) to their initial 3D
            coordinates, represented as NumPy arrays of shape (3,). May contain
            positions for a subset of vertices.
        fixed (list[Expr] or None, optional):
            List of vertices (SymPy expressions) whose positions should remain fixed
            at their coordinates specified in `pos0`. Defaults to None

    Returns:
        dict:
            A mapping from each vertex in `K` to a NumPy array of shape (3,)
            representing its computed 3D coordinates.
    """
    # Extract vertices
    vertices = [v.get_expression() for v in K[0]]

    # Assign an initial position to the vertices
    if pos0 is None:
        letters = list(K.get_alphabet().letters.values())
        pos0 = get_starting_position(letters, vertices)

    # Fixed vertices positions
    if sym.sympify(1) in vertices:
        fixed = [sym.sympify(1)] if fixed is None else fixed

    # Compute Weighted edges
    all_faces = K.get_complex(
        max_dim=5
    )  # We only use information of the 5-skeleton to plot
    edges_extra = []
    for dim in range(1, all_faces.dim + 1):
        for block in all_faces[dim]:
            new_edge = [block.min_word, block.max_word, 1/dim**3]
            edges_extra.append(new_edge)

    graph_extra = nx.Graph()
    graph_extra.add_nodes_from(vertices)
    graph_extra.add_weighted_edges_from(edges_extra)

    # Compute initial pos
    pos = pos0.copy()
    for vertex in vertices:
        if not vertex in pos:
            pos[vertex]=np.random.normal([0,0,0],[1,1,1])
    
    #pos2 = nx.kamada_kawai_layout(graph_extra, pos=pos, dim=3)
    pos2 = nx.spring_layout(graph_extra, k=1, pos=pos, dim=3, fixed=fixed)

    # Update
    pos.update(pos2)

    # If there is more than one connected component, make positions closer
    components = list(nx.connected_components(graph_extra))

    # Most external point
    component0 = components[0]
    coords = [pos[k] for k in component0]
    centroid = np.mean(coords, axis=0)
    i = np.argmax(np.sum((np.array(coords) - centroid) ** 2, axis=1))
    prev_bound = coords[i]
    for num in range(1, len(components)):
        pos, prev_bound = recenter(components[num], pos, prev_bound)

    return pos
