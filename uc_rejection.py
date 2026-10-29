"""Unit-cell clustering retained for standalone use."""
from pathlib import Path


def ucr(data, nClusters, dir_output, single=False):
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells = np.loadtxt(data, delimiter=',', ndmin=2)
    if cells.shape[0] == 0 or not np.isfinite(cells).all():
        raise ValueError('Unit cells must contain finite values')
    if nClusters < 1:
        raise ValueError('Cluster count must be positive')
    out = Path(dir_output)
    out.mkdir(parents=True, exist_ok=True)
    if len(cells) == 1:
        assignments = np.ones(1, dtype=int)
    else:
        # Retain every cell parameter, including constant crystallographic angles.
        # Standardize only dimensions with nonzero variance for single linkage.
        if single:
            scale = cells.std(axis=0)
            scale[scale == 0] = 1
            tree = linkage(cells / scale, method='single')
        else:
            tree = linkage(cells, method='ward')
        assignments = fcluster(tree, min(nClusters, len(cells)), criterion='maxclust')
        fig, ax = plt.subplots()
        dendrogram(tree, orientation='right', labels=[str(i+1) for i in range(len(cells))], ax=ax)
        fig.savefig(out / 'unit_cell_classes.pdf', bbox_inches='tight')
        plt.close(fig)
    result = [[int(c), i+1] for i, c in enumerate(assignments)]
    (out/'cluster.txt').write_text(''.join('cluster%sq %s\n' % tuple(x) for x in result))
    return result


def single(n_data, nClusters, dir_output):
    return [[nClusters, i] for i in range(1, n_data+1)]
