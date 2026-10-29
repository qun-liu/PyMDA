"""Discovery and statistics helpers independent of the CCP4 Python runtime."""
import hashlib
import re
from pathlib import Path


def discover_masters(root):
    """Skip broken links, raster scans and duplicate master-file links."""
    found = set()
    for path in Path(root).expanduser().rglob('*master.h5'):
        if path.is_file() and 'raster' not in str(path).lower():
            found.add(path.resolve())
    return sorted(found)


def dataset_id(path):
    return path.name.removesuffix('_master.h5') + '_' + hashlib.sha256(str(path).encode()).hexdigest()[:8]


def aimless_stats(path):
    text = Path(path).read_text()
    labels = {'resolution': 'High resolution limit', 'isigi': 'Mean((I)/sd(I))', 'cc_half': 'Mn(I) half-set correlation CC(1/2)', 'completeness': 'Completeness', 'multiplicity': 'Multiplicity', 'rmeas': 'Rmeas (all I+ & I-)'}
    result = {}
    for key, label in labels.items():
        matches = re.findall(r'^' + re.escape(label) + r'\s+(.+)$', text, re.M)
        if matches:
            numbers = re.findall(r'-?\d+(?:\.\d+)?', matches[-1])
            if len(numbers) >= 3:
                result[key] = [float(x) for x in numbers[:3]]
    if not all(k in result for k in ('resolution', 'isigi', 'cc_half', 'completeness')):
        raise ValueError('Missing AIMLESS summary statistics in %s' % path)
    return result


def admissible(stats, completeness=90, cc=0.3, isig=1, outer_completeness=80):
    return (stats['completeness'][0] >= completeness
            and stats['completeness'][2] >= outer_completeness
            and stats['cc_half'][2] >= cc and stats['isigi'][2] >= isig)
