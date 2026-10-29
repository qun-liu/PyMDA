"""Standalone clustering from condensed pairwise intensity correlations."""
from pathlib import Path


def pr(datafile, clusts, dir_output):
    import numpy as np
    from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    from scipy.spatial.distance import num_obs_y
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    correlations = np.fromstring(Path(datafile).read_text(), sep=' ')
    if not len(correlations) or not np.isfinite(correlations).all():
        raise ValueError('Expected finite condensed pairwise correlations')
    num_obs_y(correlations)  # Validate triangular length.
    if clusts < 1 or np.any(np.abs(correlations) > 1):
        raise ValueError('Invalid cluster count or correlation outside [-1,1]')
    tree = linkage(1-correlations, method='average')
    assignments = fcluster(tree, clusts, criterion='maxclust')
    out = Path(dir_output)
    out.mkdir(parents=True, exist_ok=True)
    result = [[int(c), i+1] for i,c in enumerate(assignments)]
    (out/'cluster.txt').write_text(''.join('clusterb%s %s\n' % tuple(x) for x in result))
    fig, ax = plt.subplots()
    dendrogram(tree, orientation='right', labels=[str(i+1) for i in range(len(assignments))], ax=ax)
    fig.savefig(out/'1-picc.pdf', bbox_inches='tight')
    plt.close(fig)
    return result
