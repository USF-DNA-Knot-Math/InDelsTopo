from scipy.cluster.hierarchy import dendrogram
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


# --- Color / label helpers --------------------------------------------------

def get_color_human(label):
    if 'AB' in label:
        return 'red'
    if 'CD' in label:
        return 'orange'
    if 'A' in label:
        return 'green'
    if 'B' in label:
        return 'blue'
    return 'black'


def get_label_human(label):
    if 'AB' in label:
        return 'AB'
    if 'CD' in label:
        return 'CD'
    if 'A' in label:
        return 'A'
    if 'B' in label:
        return 'B'
    return 'Other'


def get_color_yeast(label):
    if 'R1_antisense' in label or 'R2_antisense' in label:
        return 'red'
    if 'R1_branch' in label or 'R2_branch' in label:
        return 'blue'
    return 'black'


def get_label_yeast(label):
    if 'R1_antisense' in label or 'R2_antisense' in label:
        return 'antisense'
    if 'R1_branch' in label or 'R2_branch' in label:
        return 'branch'
    return 'Other'


# --- Plotting ---------------------------------------------------------------

def plot_dendrogram(model, color='human', **kwargs):
    """Plot dendrogram and recolor leaf labels."""
    if color == 'human':
        get_color = get_color_human
    elif color == 'yeast':
        get_color = get_color_yeast
    else:
        raise ValueError(f"Unknown color: {color}")

    # count samples under each internal node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)

    for i, merge in enumerate(model.children_):
        counts[i] = sum(
            1 if child < n_samples else counts[child - n_samples]
            for child in merge
        )

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    dendrogram(
        linkage_matrix,
        orientation='right',
        color_threshold=0,
        **kwargs
    )

    # recolor y-axis labels
    ax = plt.gca()
    for lbl in ax.get_ymajorticklabels():
        lbl.set_color(get_color(lbl.get_text()))

    plt.gcf().subplots_adjust(bottom=0.5)


def make_legend_handles(raw_names, get_label, get_color, markersize=8):
    """Create legend handles for scatter plots."""
    seen = set()
    handles = []

    for name in raw_names:
        label = get_label(name)
        if label in seen:
            continue
        seen.add(label)

        handles.append(
            Line2D(
                [0], [0],
                marker='o',
                linestyle='',
                markerfacecolor=get_color(name),
                markeredgecolor='none',
                markersize=markersize,
                label=label,
            )
        )

    return handles


def plot_projection(X, labels, color='human', figsize=(10, 10), show_text=False):
    """2D scatter projection with external legend."""
    if color == 'human':
        get_color = get_color_human
        get_label = get_label_human
    elif color == 'yeast':
        get_color = get_color_yeast
        get_label = get_label_yeast
    else:
        raise ValueError(f"Unknown color: {color}")

    plt.figure(figsize=figsize)

    for i, label in enumerate(labels):
        plt.scatter(X[i, 0], X[i, 1], color=get_color(label), s=100)
        if show_text:
            plt.text(X[i, 0], X[i, 1], label, color=get_color(label))

    handles = make_legend_handles(labels, get_label, get_color)

    plt.legend(
        handles=handles,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=1.5,
        fontsize=12
    )

# --- Utilities --------------------------------------------------------------

def safe_get(L, i, default=0):
    return L[i] if 0 <= i < len(L) else default
