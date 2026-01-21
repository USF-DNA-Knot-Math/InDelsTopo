from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np


def get_color_human(label):
    if 'AB' in label:
        return 'red'
    if 'CD' in label:
        return 'orange'
    if 'A' in label:
        return 'green'
    if 'B' in label:
        return 'blue'
    else:
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
    else:
        return 'Other'


def get_color_yeast(label):
    if 'R1_antisense' in label:
        return 'red'
    if 'R2_antisense' in label:
        return 'red'
    if 'R1_branch' in label:
        return 'blue'
    if 'R2_branch' in label:
        return 'blue'
    else:
        return 'black'

def get_label_yeast(label):
    if 'R1_antisense' in label:
        return 'antisense'
    if 'R2_antisense' in label:
        return 'antisense'
    if 'R1_branch' in label:
        return 'branch'
    if 'R2_branch' in label:
        return 'branch'
    else:
        return 'Other'

def plot_dendrogram(model, color='human', **kwargs):
    # Create linkage matrix and then plot the dendrogram

    if color == 'human':
        get_color = get_color_human
    elif color == 'yeast':
        get_color = get_color_yeast
    else:
        raise ValueError(f"Unknown color: {color}")

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, orientation='right',**kwargs,color_threshold=0)
    ax = plt.gca()
    ylbls = ax.get_ymajorticklabels()
    for i in range(len(ylbls)):
        ylbls[i].set_color(get_color(ylbls[i].get_text()))
    plt.gcf().subplots_adjust(bottom=0.5)


def make_legend_handles(raw_names, get_label, get_color, markersize=8):
    """
    Create legend handles for a scatter plot based on raw names.
    """
    seen = set()
    handles = []

    for name in raw_names:
        label = get_label(name)
        if label in seen:
            continue
        seen.add(label)

        color = get_color(name)

        handles.append(
            Line2D(
                [0], [0],
                marker='o',
                linestyle='',
                markerfacecolor=color,
                markeredgecolor='none',
                markersize=markersize,
                label=label,
            )
        )

    return handles


def plot_projection(X, labels, color='human', figsize=(10,10),show_text=False):
    if color == 'human':
        get_color = get_color_human
        get_label=get_label_human
    elif color == 'yeast':
        get_color = get_color_yeast
        get_label = get_label_yeast
    else:
        raise ValueError(f"Unknown color: {color}")
    
    plt.figure(figsize=figsize)
    for i,label in enumerate(labels):
        plt.scatter(X[i,0],X[i,1],color=get_color(label),s=100)
        if show_text:
            plt.text(X[i,0],X[i,1],label,color=get_color(label))
    #produce legend
    handles=make_legend_handles(labels, get_label, get_color)
    plt.legend(handles=handles,  loc='center left',  bbox_to_anchor=(1.02, 0.5),
               frameon=False, markerscale=1.5, fontsize=12)
    plt.show()

def safe_get(L, i, default=0):
    return L[i] if 0 <= i < len(L) else default