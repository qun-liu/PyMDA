"""Plot the measured resolution/signal/completeness tradeoff using CCP4 Python."""
import json
from pathlib import Path
import sys
import os

out = Path(sys.argv[1])
os.environ.setdefault('MPLCONFIGDIR', str(out/'plot_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
trials = json.loads((out/'trials.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(13,4), sharex=True)
for policy, color in [('selected','tab:blue'), ('cchalf','tab:green'), ('full','tab:orange')]:
    rows = [r for r in trials if r['prefix_policy']==policy]
    label = {'selected':'I/sigma prefixes', 'cchalf':'CC-half prefixes', 'full':'Full scans'}[policy]
    for ax, (key, index, title) in zip(axes, [('isigi',0,'Merged I / sigma(I), overall'), ('completeness',0,'Completeness, overall (%)'), ('completeness',2,'Completeness, outer shell (%)')]):
        ax.scatter([r['statistics']['resolution'][0] for r in rows], [r['statistics'][key][index] for r in rows], s=18, alpha=.5, color=color, label=label)
        ax.set_title(title)
        ax.set_xlabel('High-resolution cutoff (A)')
        ax.grid(alpha=.2)
        ax.invert_xaxis() if not ax.xaxis_inverted() else None
best = json.loads((out/'best.json').read_text())
for ax, (key,index) in zip(axes,[('isigi',0),('completeness',0),('completeness',2)]):
    ax.scatter([best['statistics']['resolution'][0]],[best['statistics'][key][index]],marker='*',s=180,color='black',label='chosen')
axes[1].axhline(best['thresholds']['completeness_overall'], color='gray', linestyle='--', linewidth=1)
axes[2].axhline(best['thresholds']['completeness_outer'], color='gray', linestyle='--', linewidth=1)
axes[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out/'resolution_tradeoff.png',dpi=180)
plt.close(fig)
