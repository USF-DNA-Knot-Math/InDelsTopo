import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE
from scipy.cluster.hierarchy import dendrogram

# Color helper
# --------------------------
def get_color(file_name):
    """Assign colors based on experiment name."""
    if 'LIN' in file_name:
        return 'red'
    elif 'SUP' in file_name:
        return 'blue'
    else:
        return 'green'

# --------------------------
# Plot Euler curves
# --------------------------
def plot_euler_curves(exp_names, euler_curves, title=None):
    fig, axes = plt.subplots(3, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, ax in enumerate(axes):
        x, y = euler_curves[i]
        ax.plot(x, y)
        ax.set_title(exp_names[i])

    if title:
        fig.suptitle(title, fontsize=16)  # Add group title above all subplots

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for suptitle
    plt.show()

# --------------------------
# t-SNE visualization
# --------------------------
def plot_tsne(exp_names, euler_curves, perplexity=2, random_state=13, metric='cosine', figsize=(5,5), title=None):
    X = np.array([curve[1] for curve in euler_curves])
    tSNE = TSNE(metric=metric, random_state=random_state, perplexity=perplexity, n_components=2)
    Y = tSNE.fit_transform(X)

    plt.figure(figsize=figsize)
    plt.scatter(Y[:,0], Y[:,1])

    for i, label in enumerate(exp_names):
        plt.text(Y[i,0], Y[i,1], label, color=get_color(label))
    
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    if title:
        plt.title(title)
    plt.show()



def plot_dendrogram(model, labels=None, **kwargs):
    """
    Plot a dendrogram from a fitted AgglomerativeClustering model.
    """
    n_samples = len(model.labels_)

    # Count samples under each internal node
    counts = np.zeros(model.children_.shape[0])
    for i, merge in enumerate(model.children_):
        counts[i] = sum(
            1 if child < n_samples else counts[child - n_samples]
            for child in merge
        )

    # Create linkage matrix for dendrogram
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    # If labels are not provided, default to indices
    if labels is None:
        labels = [str(i) for i in range(n_samples)]

    # Plot dendrogram
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        orientation='right',
        color_threshold=0,
        **kwargs
    )

    # Recolor y-axis labels if function provided
    ax = plt.gca()
    for lbl in ax.get_ymajorticklabels():
        text = lbl.get_text()
        lbl.set_color(get_color(text))

    # Adjust subplot to avoid clipping
    plt.gcf().subplots_adjust(bottom=0.5)
