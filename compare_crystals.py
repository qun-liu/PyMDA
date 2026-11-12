"""Run with CCP4 Python to compare cells and aligned common intensities."""
import json
import math
from pathlib import Path
import shutil
import sys
from iotbx import mtz
from autoDials import Toolchain


def intensity(path, cutoff):
    arrays = mtz.object(file_name=str(path)).as_miller_arrays()
    candidates = [a for a in arrays if a.is_xray_intensity_array() and not a.anomalous_flag()]
    if not candidates:
        candidates = [a.average_bijvoet_mates() for a in arrays if a.is_xray_intensity_array()]
    if not candidates:
        raise ValueError('No intensity array in MTZ')
    array = candidates[0].resolution_filter(d_min=cutoff)
    return array.merge_equivalents().array().map_to_asu()


def compare(source, ccp4, cutoff, length_tolerance, angle_tolerance, min_cc, min_common, dataset_ids=None):
    tools = Toolchain(ccp4)
    records = [(p.parent.name, json.loads(p.read_text())) for p in sorted(source.glob('crystals/*/manifest.json'))]
    if dataset_ids:
        records = [(key,record) for key,record in records if key in dataset_ids]
    datasets = []
    for identifier, record in records:
        array = intensity(Path(record['selected']['mtz']), cutoff)
        reduced = array.crystal_symmetry().niggli_cell().unit_cell().parameters()
        datasets.append({'id': identifier, 'name': Path(record['config']['master']).name,
                         'cell': list(array.unit_cell().parameters()), 'reduced_cell': list(reduced),
                         'last_frame': record['selected'].get('last_frame'), 'mtz': str(Path(record['selected']['mtz']).resolve()), 'space_group': str(array.space_group_info())})
    pairs = []
    for i in range(len(records)):
        for j in range(i+1, len(records)):
            a,b = datasets[i],datasets[j]
            ca,cb = a['reduced_cell'],b['reduced_cell']
            lengths = max(abs(x-y)/((x+y)/2)*100 for x,y in zip(ca[:3],cb[:3]))
            angles = max(abs(x-y) for x,y in zip(ca[3:],cb[3:]))
            pair = {'a': a['id'], 'b': b['id'], 'length_difference_percent': lengths,
                    'angle_difference_degrees': angles, 'cell_compatible': lengths <= length_tolerance and angles <= angle_tolerance,
                    'correlation': None, 'common_reflections': 0, 'error': ''}
            try:
                ref = Path(records[i][1]['selected']['mtz'])
                target = Path(records[j][1]['selected']['mtz'])
                work = source/'comparisons'/f'{i+1}_{j+1}'
                work.mkdir(parents=True, exist_ok=True)
                # Pointless splits command-line filenames containing spaces.
                # Local names preserve access to exports in any source folder.
                for name, original in [('reference.mtz',ref),('target.mtz',target)]:
                    local = work/name
                    local.unlink(missing_ok=True)
                    try:
                        local.symlink_to(original.resolve())
                    except OSError:
                        shutil.copy2(original,local)
                tools.run('pointless', ['hklin', 'target.mtz', 'hklref', 'reference.mtz', 'hklout', 'aligned.mtz'], work, 'END\n', ['aligned.mtz'])
                x,y = intensity(ref,cutoff),intensity(work/'aligned.mtz',cutoff)
                x,y = x.common_sets(y, assert_is_similar_symmetry=True)
                pair['common_reflections'] = x.size()
                if x.size() < min_common:
                    raise ValueError('Insufficient common reflections')
                cc = x.correlation(y, assert_is_similar_symmetry=True).coefficient()
                if cc is None or not math.isfinite(cc):
                    raise ValueError('Correlation undefined (constant intensities)')
                pair['correlation'] = cc
            except Exception as error:
                pair['error'] = str(error)
            pair['compatible'] = pair['cell_compatible'] and pair['correlation'] is not None and pair['correlation'] >= min_cc
            pairs.append(pair)
    # Complete-link groups require every pair to pass; no transitive bridging.
    groups = []
    lookup = {frozenset((p['a'],p['b'])):p['compatible'] for p in pairs}
    for d in datasets:
        for group in groups:
            if all(lookup.get(frozenset((d['id'],other)),False) for other in group):
                group.append(d['id']); break
        else:
            groups.append([d['id']])
    result = {'datasets':datasets,'pairs':pairs,'groups':groups,'criteria':{'resolution':cutoff,'length_percent':length_tolerance,'angle_degrees':angle_tolerance,'min_correlation':min_cc,'min_common':min_common},'method':'Niggli-reduced cell differences; Pointless reference indexing; Pearson correlation of merged common non-anomalous intensities; deterministic complete-link groups'}
    temporary = source/'comparison.tmp'
    temporary.write_text(json.dumps(result,indent=2))
    temporary.replace(source/'comparison.json')
    return result

if __name__ == '__main__':
    compare(Path(sys.argv[1]),sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]),int(sys.argv[7]),sys.argv[8:])
