#!/usr/bin/env python3
"""Microcrystal data assembly with current DIALS and CCP4.

Compare cumulative frame prefixes and crystal subsets; retain all trials and
select supported resolution first, then merged I/sigma(I).
"""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
from autoDials import Toolchain, default_ccp4
from helpers import discover_masters, dataset_id, frame_wedges, aimless_stats, admissible


def options():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset_ids', nargs='+', help='Choose processed dataset IDs for comparison or assembly')
    p.add_argument('--compare', action='store_true')
    p.add_argument('--cell_length_tolerance', type=float, default=2)
    p.add_argument('--cell_angle_tolerance', type=float, default=2)
    p.add_argument('--min_pair_cc', type=float, default=0.9)
    p.add_argument('--min_common', type=int, default=100)
    detected = default_ccp4()
    p.add_argument('--ccp4', default=str(detected) if detected else None,
                   help='CCP4 installation; defaults to PYMDA_CCP4 or the installation found by which aimless')
    p.add_argument('--hdf5_path', default='~/.codex/pymda_data')
    p.add_argument('--master_files', nargs='+', help='Process only these discovered rotation master files')
    p.add_argument('--dir_output', default='~/.codex/pymda_output')
    p.add_argument('--run_dials', action='store_true')
    p.add_argument('--run_mda', action='store_true')
    p.add_argument('--discover', action='store_true')
    p.add_argument('--exports_only', action='store_true', help='Recover exports and merging from completed integrations without reading raw images')
    p.add_argument('--thread', type=int, default=2)
    p.add_argument('--merge_workers', type=int, help='Independent merging jobs; defaults to thread count, capped at four')
    p.add_argument('--subset_search', choices=['auto', 'exhaustive', 'bounded'], default='auto', help='Auto compares all subsets through ten datasets and uses a bounded search for larger collections')
    p.add_argument('--max_subsets', type=int, default=64, help='Maximum candidate subsets in a bounded search')
    p.add_argument('--wedges', type=int, default=4, help='Number of cumulative frame prefixes per crystal')
    p.add_argument('--split_wedges', type=int, default=1, help='Consecutive wedges per scan; 3 splits 1800 frames into 1-600, 601-1200, 1201-1800')
    p.add_argument('--frames', type=int, help='Optional maximum frames per crystal')
    p.add_argument('--reso', type=float, default=1.8, help='Integration high-resolution limit, A')
    p.add_argument('--reso_cchalf', type=float, default=2.5, help='Common resolution for prefix selection')
    p.add_argument('--resolutions', type=float, nargs='+', default=[1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.0])
    p.add_argument('--spg', help='Known space group; otherwise index without imposing symmetry')
    p.add_argument('--unit_cell')
    p.add_argument('--fix_detector', action='store_true')
    p.add_argument('--min_completeness', type=float, default=90, help='Overall completeness threshold (paper uses over 90 percent)')
    p.add_argument('--min_outer_completeness', type=float, default=80)
    p.add_argument('--prefix_policy', choices=['selected', 'cchalf', 'full', 'both', 'all'], default='all', help='Compare I/sigma-selected, CC-half-selected and full-scan prefixes')
    p.add_argument('--min_cc_half', type=float, default=0.3)
    p.add_argument('--min_isigi', type=float, default=1.0, help='Minimum outer-shell merged I/sigma(I)')
    p.add_argument('--nres', type=int, help='Known number of residues in the asymmetric unit')
    p.add_argument('--crystal_identity', choices=['unknown', 'repeat', 'independent'], default='unknown', help='Identity established by collection metadata; not inferred from unit cells')
    return p.parse_args()


def save(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2))
    temp.replace(path)


def write_report(out, best, alternatives, records, identity='unknown'):
    lines = ['# PyMDA processing results', '',
             'Processed %s rotation datasets with current DIALS and CCP4.' % len(records), '',
             '| Selection | Resolution (A) | I/sigma overall | I/sigma outer | CC-half outer | Completeness overall / outer (%) | Datasets |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for label, result in [('Highest supported resolution', best)] + list(alternatives.items()):
        s = result['statistics']
        lines.append('| %s | %.2f | %.1f | %.1f | %.3f | %.1f / %.1f | %s |' % (label, s['resolution'][0], s['isigi'][0], s['isigi'][2], s['cc_half'][2], s['completeness'][0], s['completeness'][2], len(result['crystals'])))
    lines += ['', 'Final amplitudes and CCP4 Free-R flags: `best/final.mtz`.',
              'Merged intensities and anomalous columns: `best/scaled.mtz`.',
              'Thresholds: `%s`.' % json.dumps(best['thresholds'], sort_keys=True), '',
              'Chosen input scans and last included frame:', '']
    lines += ['- `%s`: frame %s' % pair for pair in zip(best['crystals'], best['last_frames'])]
    identity_note = {'unknown': 'Crystal identity is unverified. This output establishes a multi-dataset merge and does not establish that the inputs are independent crystals.', 'repeat': 'The inputs were identified as repeat scans of one crystal. This output is a multi-dataset merge.', 'independent': 'The inputs were identified as independent crystals by the supplied collection information.'}[identity]
    if len(records) == 1:
        identity_note = 'This result uses one rotation dataset.'
    labels = {Path(r['config']['master']).parent.name for r in records}
    if identity == 'unknown' and len(labels) == 1:
        identity_note += ' All input paths use the same immediate folder label: `%s`.' % next(iter(labels))
    lines += ['', 'The search compares cumulative prefixes chosen by I/sigma, prefixes chosen by CC-half, and full scans, across the selected dataset subsets and the recorded cutoff grid. It does not exhaust all possible frame combinations or reproduce the paper\'s SmRmerge rejection procedure.', '', identity_note, '',
              'See `trials.csv`, `trials.json`, `resolution_tradeoff.png`, and individual tool logs for the full comparison. CCP4 Free-R value 0 identifies the approximately 5% test set; other values identify working subsets.']
    if (out/'selection.json').is_file():
        lines += ['', 'Dataset selection provenance: `selection.json`.']
    if (out/'comparison.json').is_file():
        lines += ['', 'Unit-cell and common-intensity screening: `comparison.json`.']
    if (out/'search.json').is_file():
        search = json.loads((out/'search.json').read_text())
        lines += ['', 'Subset search: %s; %s candidates for %s datasets. See `search.json` for the exact subsets and method.' %
                  (search['mode'], search['candidate_count'], search['dataset_count'])]
    for record in records:
        if record.get('scan_wedges'):
            lines += ['', 'Consecutive wedges for `%s`: %s. Individual MTZs, actual integrated bounds and statistics are recorded in its crystal manifest.' %
                      (record['config']['master'], ', '.join('%s-%s' % (w['first_frame'], w['last_frame']) for w in record['scan_wedges']))]
    (out/'processing_report.md').write_text('\n'.join(lines)+'\n')


def scale(tools, mtzs, directory, resolution):
    directory.mkdir(parents=True, exist_ok=True)
    summary = directory / 'statistics.json'
    inputs = {'mtzs': [str(x.resolve()) for x in mtzs], 'resolution': resolution}
    newest_input = max(p.stat().st_mtime_ns for p in mtzs)
    if summary.is_file() and summary.stat().st_mtime_ns >= newest_input and (directory / 'scaled.mtz').is_file() and (directory / 'inputs.json').is_file() and json.loads((directory / 'inputs.json').read_text()) == inputs:
        return json.loads(summary.read_text())
    if summary.exists():
        summary.unlink()
    save(directory / 'inputs.json', inputs)
    commands = '\n'.join('HKLIN "%s"' % x.resolve() for x in mtzs)
    # Let Pointless check symmetry and align indexing; do not blindly COPY.
    reusable = None
    for sibling in sorted(directory.parent.glob('resolution_*')):
        metadata = sibling / 'inputs.json'
        if sibling != directory and metadata.is_file() and (sibling/'sorted.mtz').is_file() and (sibling/'statistics.json').is_file() and (sibling/'statistics.json').stat().st_mtime_ns >= newest_input and json.loads(metadata.read_text())['mtzs'] == inputs['mtzs']:
            reusable = sibling
            break
    if reusable:
        for name in ['sorted.mtz', 'pointless.console.log', 'pointless.command.json']:
            shutil.copy2(reusable/name, directory/name)
        save(directory/'pointless_reused.json', {'source': str(reusable)})
    else:
        tools.run('pointless', ['hklout', 'sorted.mtz'], directory,
                  commands + '\nALLOW OUTOFSEQUENCEFILES\nEND\n', ['sorted.mtz'])
    log = tools.run('aimless', ['hklin', 'sorted.mtz', 'hklout', 'scaled.mtz'], directory,
                   'BINS 20\nRESOLUTION HIGH %.4f\nANOMALOUS ON\nEND\n' % resolution, ['scaled.mtz'])
    stats = aimless_stats(log)
    save(summary, stats)
    return stats


def export_scan_wedges(tools, directory, args):
    """Export disjoint ranges from one integration, retaining requested bounds."""
    count = getattr(args, 'split_wedges', 1)
    if count == 1:
        return []
    imported = json.loads((directory/'imported.expt').read_text())
    first, last = imported['scan'][0]['image_range']
    if args.frames:
        last = min(last, first + args.frames - 1)
    ranges = frame_wedges(first, last, count)
    integrated = json.loads((directory/'integrated.expt').read_text())
    actual_first, actual_last = integrated['scan'][0]['image_range']
    exports = []
    for number, (first, last) in enumerate(ranges, 1):
        used_first, used_last = max(first, actual_first), min(last, actual_last)
        wedge = directory/('wedge_%04d_%04d' % (first, last))
        item = {'number': number, 'first_frame': first, 'last_frame': last,
                'used_first_frame': used_first if used_first <= used_last else None,
                'used_last_frame': used_last if used_first <= used_last else None,
                'mtz': None, 'statistics': None, 'error': ''}
        try:
            if used_first > used_last:
                raise ValueError('No integrated frames remain in this wedge')
            tools.run('dials.slice_sequence', [directory/'integrated.expt', directory/'integrated.refl',
                      'image_range=%s,%s' % (used_first, used_last), 'output.experiments=sliced.expt',
                      'output.reflections=sliced.refl'], wedge, expected=['sliced.expt', 'sliced.refl'])
            tools.run('dials.export', ['sliced.expt', 'sliced.refl', 'mtz.hklout=integrated.mtz'],
                      wedge, expected=['integrated.mtz'])
            stats = scale(tools, [wedge/'integrated.mtz'], wedge/'assessment', args.reso_cchalf)
            item.update(mtz=str((wedge/'integrated.mtz').resolve()), statistics=stats)
            if (wedge/'failure.json').exists():
                (wedge/'failure.json').unlink()
        except (RuntimeError, ValueError) as error:
            item['error'] = str(error)
            wedge.mkdir(parents=True, exist_ok=True)
            save(wedge/'failure.json', item)
        exports.append(item)
    return exports


def record_scan_wedges(tools, directory, record, args):
    record['split_wedges'] = getattr(args, 'split_wedges', 1)
    record['scan_wedges'] = export_scan_wedges(tools, directory, args)
    save(directory/'manifest.json', record)
    if record['scan_wedges'] and all(w['error'] for w in record['scan_wedges']):
        raise RuntimeError('No consecutive wedge processed successfully; inspect %s' % (directory/'manifest.json'))
    return record


def process(tools, master, out, args):
    directory = out / 'crystals' / dataset_id(master)
    manifest = directory / 'manifest.json'
    config = {'master': str(master), 'frames': args.frames, 'wedges': args.wedges, 'resolution': args.reso, 'space_group': args.spg, 'unit_cell': args.unit_cell, 'fix_detector': args.fix_detector}
    if manifest.is_file():
        old = json.loads(manifest.read_text())
        if old['config'] != config:
            raise ValueError('Processing settings changed; use a new output directory: %s' % directory)
        if all(Path(p['mtz']).is_file() for p in old['prefixes']):
            for prefix in old['prefixes']:
                mtz = Path(prefix['mtz'])
                prefix['statistics'] = scale(tools, [mtz], mtz.parent/'assessment', args.reso_cchalf)
            old['selected'] = max(old['prefixes'], key=lambda x: (x['statistics']['cc_half'][0] >= 0.3, x['statistics']['isigi'][0], x['statistics']['cc_half'][0]))
            return record_scan_wedges(tools, directory, old, args)
    directory.mkdir(parents=True, exist_ok=True)
    config_file = directory / 'processing_config.json'
    had_config = config_file.is_file()
    if had_config and json.loads(config_file.read_text()) != config:
        raise ValueError('Processing settings changed; use a new output directory: %s' % directory)
    save(config_file, config)
    # A completed integration can be recovered after an export/scaling failure.
    # Verify its command, source import and completion log before adopting it.
    integration_command = directory / 'dials.integrate.command.json'
    import_command = directory / 'dials.import.command.json'
    integrated_log = directory / 'dials.integrate.console.log'
    recover = all((directory / f).is_file() for f in ['integrated.expt', 'integrated.refl', 'imported.expt']) and integration_command.is_file() and import_command.is_file() and integrated_log.is_file()
    if recover:
        argv = json.loads(integration_command.read_text())['argv']
        imported_argv = json.loads(import_command.read_text())['argv']
        first, last = json.loads((directory/'imported.expt').read_text())['scan'][0]['image_range']
        if args.frames:
            last = min(last, first + args.frames - 1)
        recover = (str(master) in imported_argv and 'prediction.d_min=%s' % args.reso in argv
                   and 'Saving the experiments to integrated.expt' in integrated_log.read_text()
                   and 'scan_range=%s,%s' % (first, last) in argv
                   and (had_config or (args.spg is None and args.unit_cell is None and not args.fix_detector)))
    if args.exports_only and not recover:
        raise ValueError('Cannot safely recover this integration; run normal processing: %s' % directory)
    if not recover:
        tools.run('dials.import', [master], directory, expected=['imported.expt'])
    exp = json.loads((directory / 'imported.expt').read_text())
    start, end = exp['scan'][0]['image_range']
    if args.frames:
        end = min(end, start + args.frames - 1)
    frame_wedges(start, end, getattr(args, 'split_wedges', 1))
    if not recover:
        tools.run('dials.find_spots', ['imported.expt', 'scan_range=%s,%s' % (start, end), 'filter.d_min=%s' % args.reso, 'nproc=%s' % args.thread], directory, expected=['strong.refl'])
    index_args = ['imported.expt', 'strong.refl']
    if args.spg:
        index_args.append('space_group=%s' % args.spg)
    if args.unit_cell:
        index_args.append('unit_cell=%s' % args.unit_cell)
    if args.fix_detector:
        index_args.append('detector.fix=all')
    if not recover:
        tools.run('dials.index', index_args, directory, expected=['indexed.expt', 'indexed.refl'])
        refine_args = ['indexed.expt', 'indexed.refl']
        if args.fix_detector:
            refine_args.append('detector.fix=all')
        tools.run('dials.refine', refine_args, directory, expected=['refined.expt', 'refined.refl'])
        tools.run('dials.integrate', ['refined.expt', 'refined.refl', 'prediction.d_min=%s' % args.reso, 'scan_range=%s,%s' % (start, end), 'nproc=%s' % args.thread], directory, expected=['integrated.expt', 'integrated.refl'])
    # DIALS may trim scan endpoints lacking reflections during refinement.
    integrated = json.loads((directory / 'integrated.expt').read_text())
    actual_start, actual_end = integrated['scan'][0]['image_range']
    start, end = max(start, actual_start), min(end, actual_end)
    prefixes = []
    for last in sorted(set(start - 1 + ((end-start+1)*i + args.wedges-1)//args.wedges for i in range(1, args.wedges+1))):
        wedge = directory / ('prefix_%04d' % last)
        tools.run('dials.slice_sequence', [directory / 'integrated.expt', directory / 'integrated.refl', 'image_range=%s,%s' % (start,last), 'output.experiments=sliced.expt', 'output.reflections=sliced.refl'], wedge, expected=['sliced.expt', 'sliced.refl'])
        tools.run('dials.export', ['sliced.expt', 'sliced.refl', 'mtz.hklout=integrated.mtz'], wedge, expected=['integrated.mtz'])
        try:
            stats = scale(tools, [wedge / 'integrated.mtz'], wedge / 'assessment', args.reso_cchalf)
            prefixes.append({'last_frame': last, 'mtz': str((wedge/'integrated.mtz').resolve()), 'statistics': stats})
        except (RuntimeError, ValueError) as e:
            save(wedge / 'failure.json', {'error': str(e)})
    if not prefixes:
        raise RuntimeError('No prefix scaled successfully: %s' % master)
    # At a common cutoff, choose strongest supported signal, allowing partial completeness.
    best = max(prefixes, key=lambda x: (x['statistics']['cc_half'][0] >= 0.3, x['statistics']['isigi'][0], x['statistics']['cc_half'][0]))
    value = {'config': config, 'prefixes': prefixes, 'selected': best}
    return record_scan_wedges(tools, directory, value, args)


def record_id(record):
    return record.get('dataset_id') or dataset_id(Path(record['config']['master']))


def merge_subsets(records, out, args):
    """Choose an auditable finite search for a collection of any size."""
    count = len(records)
    if not count:
        raise ValueError('Processing requires at least one successful dataset')
    mode = getattr(args, 'subset_search', 'auto')
    limit = getattr(args, 'max_subsets', 64)
    if mode == 'exhaustive' and count > 10:
        raise ValueError('Exhaustive search supports at most ten datasets; choose auto or bounded')
    if count == 1:
        mode, subsets = 'single', [(0,)]
    elif mode == 'exhaustive' or (mode == 'auto' and count <= 10):
        mode = 'exhaustive'
        subsets = [indices for size in range(2, count+1)
                   for indices in itertools.combinations(range(count), size)]
    else:
        mode, subsets, seen = 'bounded', [], set()
        def add(indices):
            candidate = tuple(sorted(indices))
            if len(candidate) >= 2 and candidate not in seen and len(subsets) < limit:
                subsets.append(candidate)
                seen.add(candidate)
        add(range(count))
        comparison = out/'comparison.json'
        if comparison.is_file():
            groups = json.loads(comparison.read_text()).get('groups', [])
            identifiers = {record_id(r): i for i,r in enumerate(records)}
            for group in sorted(groups, key=lambda g: (-len(g), g)):
                add(identifiers[key] for key in group if key in identifiers)
        def quality(index):
            stats = records[index]['selected']['statistics']
            return (-(stats['cc_half'][0] >= 0.3), -stats['isigi'][0],
                    -stats['cc_half'][0], -stats['completeness'][0], records[index]['config']['master'])
        ranked = sorted(range(count), key=quality)
        for excluded in reversed(ranked[-16:]):
            add(i for i in range(count) if i != excluded)
        for size in range(2, min(count, 32)+1):
            add(ranked[:size])
        size = 64
        while size < count:
            add(ranked[:size])
            size *= 2
    description = ('All dataset subsets of size two or larger.' if mode == 'exhaustive' else
                   'One dataset.' if mode == 'single' else
                   'Full collection, available compatibility groups, removal of up to sixteen individually weakest datasets, and growing subsets ranked by individual I/sigma with CC-half support; candidate budget applies in that order.')
    search = {'mode': mode, 'dataset_count': count, 'candidate_count': len(subsets),
              'max_subsets': limit if mode == 'bounded' else None, 'method': description,
              'dataset_subsets': [[record_id(records[i]) for i in indices] for indices in subsets],
              'subsets': [[records[i]['config']['master'] for i in indices] for indices in subsets]}
    save(out/'search.json', search)
    return subsets, search


def subset_label(records, indices):
    if len(records) <= 10:
        return 'crystals_' + '_'.join(str(i+1) for i in indices)
    masters = sorted(records[i]['config']['master'] for i in indices)
    if len({r['config']['master'] for r in records}) != len(records):
        masters = sorted((record_id(records[i]), records[i]['config']['master']) for i in indices)
    digest = hashlib.sha256(json.dumps(masters).encode()).hexdigest()[:16]
    return 'datasets_%s_%s' % (len(indices), digest)


def assemble(tools, records, out, args):
    subsets, search = merge_subsets(records, out, args)
    results = []
    policies = ['selected', 'cchalf', 'full'] if args.prefix_policy == 'all' else ['selected', 'full'] if args.prefix_policy == 'both' else [args.prefix_policy]
    evaluated = {}
    def evaluate_subset(policy, indices):
        rows = []
        if policy == 'selected':
            prefixes = [records[i]['selected'] for i in indices]
        elif policy == 'cchalf':
            prefixes = [max(records[i]['prefixes'], key=lambda p: (p['statistics']['cc_half'][0], p['statistics']['cc_half'][2], p['statistics']['isigi'][0])) for i in indices]
        else:
            prefixes = [max(records[i]['prefixes'], key=lambda p: p['last_frame']) for i in indices]
        mtzs = [Path(p['mtz']) for p in prefixes]
        for resolution in sorted(set(args.resolutions)):
            if resolution < args.reso:
                continue
            directory = out / 'trials' / policy / subset_label(records, indices) / ('resolution_%.2f' % resolution)
            try:
                key = (tuple(map(str, mtzs)), resolution)
                if key in evaluated:
                    stats, directory = evaluated[key]
                else:
                    stats = scale(tools, mtzs, directory, resolution)
                    evaluated[key] = (stats, directory)
                rows.append({'crystals': [records[i]['config']['master'] for i in indices], 'dataset_ids': [record_id(records[i]) for i in indices], 'prefix_policy': policy, 'last_frames': [p['last_frame'] for p in prefixes], 'requested_resolution': resolution, 'directory': str(directory), 'statistics': stats, 'eligible': admissible(stats, args.min_completeness, args.min_cc_half, args.min_isigi, args.min_outer_completeness)})
            except (RuntimeError, ValueError) as e:
                save(directory / 'failure.json', {'error': str(e)})
        return rows

    # Subsets are independent; keep resolution trials within each subset serial
    # so the Pointless output can be reused safely. Policies form separate waves.
    merge_workers = getattr(args, 'merge_workers', None) or args.thread
    with ThreadPoolExecutor(max_workers=min(merge_workers, 4)) as executor:
        for policy in policies:
            futures = [executor.submit(evaluate_subset, policy, indices) for indices in subsets]
            for future in as_completed(futures):
                results.extend(future.result())
                save(out/'trials.json', results)
    results.sort(key=lambda r: (policies.index(r['prefix_policy']), r['crystals'], r['requested_resolution']))
    save(out / 'trials.json', results)
    with (out/'trials.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['prefix_policy', 'crystals', 'last_frames', 'requested_resolution_A', 'measured_resolution_A', 'isigi_overall', 'isigi_outer', 'cc_half_outer', 'completeness_overall', 'completeness_outer', 'eligible'])
        for r in results:
            s = r['statistics']
            writer.writerow([r['prefix_policy'], ';'.join(r['crystals']), ';'.join(map(str,r['last_frames'])), r['requested_resolution'], s['resolution'][0], s['isigi'][0], s['isigi'][2], s['cc_half'][2], s['completeness'][0], s['completeness'][2], r['eligible']])
    eligible = [x for x in results if x['eligible']]
    if not eligible:
        raise RuntimeError('No trial meets the stated quality thresholds; inspect trials.json')
    best = min(eligible, key=lambda x: (x['statistics']['resolution'][0], -x['statistics']['isigi'][0], -x['statistics']['cc_half'][2]))
    best['subset_search'] = {key: search[key] for key in ('mode','dataset_count','candidate_count','max_subsets')}
    best['thresholds'] = {'completeness_overall': args.min_completeness, 'completeness_outer': args.min_outer_completeness, 'cc_half_outer': args.min_cc_half, 'isigi_outer': args.min_isigi}
    alternatives = {'strongest_signal': max(eligible, key=lambda x: (x['statistics']['isigi'][0], -x['statistics']['resolution'][0]))}
    extended = [x for x in results if x['statistics']['cc_half'][2] >= args.min_cc_half and x['statistics']['isigi'][2] >= args.min_isigi]
    if extended:
        alternatives['extended_resolution'] = min(extended, key=lambda x: (x['statistics']['resolution'][0], -x['statistics']['isigi'][0]))
    pareto = []
    strongest = float('-inf')
    for candidate in sorted(eligible, key=lambda x: (x['statistics']['resolution'][0], -x['statistics']['isigi'][0])):
        if candidate['statistics']['isigi'][0] > strongest:
            pareto.append(candidate)
            strongest = candidate['statistics']['isigi'][0]
    save(out/'pareto.json', pareto)
    for candidate in pareto:
        target = out/'alternatives'/'pareto'/('resolution_%.2f' % candidate['statistics']['resolution'][0])
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(candidate['directory'])/'scaled.mtz', target/'scaled.mtz')
        save(target/'selection.json', candidate)
    complete = [x for x in eligible if x['statistics']['completeness'][0] >= 95]
    if complete:
        alternatives['completeness_95'] = min(complete, key=lambda x: (x['statistics']['resolution'][0], -x['statistics']['isigi'][0]))
    for label, alternative in alternatives.items():
        target = out/'alternatives'/label
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(alternative['directory'])/'scaled.mtz', target/'scaled.mtz')
        save(target/'selection.json', alternative)
    save(out/'alternatives.json', alternatives)
    best['crystal_identity'] = getattr(args, 'crystal_identity', 'unknown')
    write_report(out, best, alternatives, records, best['crystal_identity'])
    final = out / 'best'
    final.mkdir(exist_ok=True)
    provenance = {'scaled_sha256': hashlib.sha256((Path(best['directory'])/'scaled.mtz').read_bytes()).hexdigest(), 'nres': args.nres, 'final_format': 'mean_and_anomalous_v1'}
    previous = final/'final_inputs.json'
    if previous.is_file() and json.loads(previous.read_text()) == provenance and (final/'final.mtz').is_file():
        save(out/'best.json', best)
        tools.run('ccp4-python', [Path(__file__).resolve().parent/'plot_results.py', out, out/'trials.json', out/'best.json'], out, expected=['resolution_tradeoff.png'])
        print(json.dumps(best, indent=2), flush=True)
        return
    for name in ('scaled.mtz', 'sorted.mtz', 'aimless.console.log', 'pointless.console.log', 'statistics.json', 'inputs.json'):
        shutil.copy2(Path(best['directory']) / name, final / name)
    trunc_args = ['-mtzin', 'scaled.mtz', '-mtzout', 'amplitudes.mtz', '-colin', '/*/*/[IMEAN,SIGIMEAN]', '-colano', '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]']
    if args.nres:
        trunc_args += ['-nres', args.nres]
    tools.run('ctruncate', trunc_args, final, expected=['amplitudes.mtz'])
    tools.run('freerflag', ['hklin', 'amplitudes.mtz', 'hklout', 'final.mtz'], final, 'FREERFRAC 0.05\nEND\n', ['final.mtz'])
    save(previous, provenance)
    save(out / 'best.json', best)
    tools.run('ccp4-python', [Path(__file__).resolve().parent/'plot_results.py', out, out/'trials.json', out/'best.json'], out, expected=['resolution_tradeoff.png'])
    print(json.dumps(best, indent=2), flush=True)


def main():
    args = options()
    if not (math.isfinite(args.cell_length_tolerance) and math.isfinite(args.cell_angle_tolerance) and args.cell_length_tolerance > 0 and args.cell_angle_tolerance > 0 and -1 <= args.min_pair_cc <= 1 and args.min_common >= 3):
        raise ValueError('Invalid crystal comparison criteria')
    if args.thread < 1 or args.wedges < 1 or (args.frames is not None and args.frames < 1) or (args.merge_workers is not None and args.merge_workers < 1):
        raise ValueError('Thread, wedge and frame counts must be positive')
    if not 1 <= args.split_wedges <= 100:
        raise ValueError('Consecutive wedge count must be between one and 100')
    if args.max_subsets < 1:
        raise ValueError('Bounded search candidate count must be positive')
    masters = discover_masters(args.hdf5_path)
    if args.master_files:
        requested = {Path(p).expanduser().resolve() for p in args.master_files}
        if not requested.issubset(set(masters)):
            raise ValueError('Selected master files must be discovered below --hdf5_path')
        masters = sorted(requested)
    if args.discover:
        print('\n'.join(map(str, masters)))
        return
    if not args.run_dials and not args.run_mda and not args.compare:
        raise ValueError('Choose --run_dials and/or --run_mda, or --discover')
    tools = Toolchain(args.ccp4)
    args.ccp4 = str(tools.root)
    out = Path(args.dir_output).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    source_dir = Path(__file__).resolve().parent
    hashes = {name: hashlib.sha256((source_dir/name).read_bytes()).hexdigest() for name in ['pymda.py', 'autoDials.py', 'helpers.py', 'plot_results.py', 'compare_crystals.py']}
    save(out / 'run.json', {'started_utc': datetime.now(timezone.utc).isoformat(), 'options': vars(args), 'ccp4': str(tools.root), 'masters': list(map(str, masters)), 'source_sha256': hashes})
    records = []
    failures = []
    if args.run_dials:
        if not masters:
            raise ValueError('No rotation master files found')
        for master in masters:
            if args.exports_only and not (out/'crystals'/dataset_id(master)/'integrated.expt').is_file():
                continue
            try:
                record = process(tools, master, out, args)
                records.append(record)
                failures.extend({'master': str(master), 'frames': [w['first_frame'], w['last_frame']],
                                 'error': w['error']} for w in record.get('scan_wedges', []) if w['error'])
            except (RuntimeError, ValueError) as e:
                print(str(e), file=sys.stderr, flush=True)
                failures.append({'master': str(master), 'error': str(e)})
                save(out / 'failures.json', failures)
    else:
        for path in sorted((out/'crystals').glob('*/manifest.json')):
            record = json.loads(path.read_text())
            record['dataset_id'] = path.parent.name
            records.append(record)
    if args.dataset_ids:
        available = {record_id(r): r for r in records}
        if not set(args.dataset_ids).issubset(available):
            raise ValueError('Selected dataset IDs must have completed manifests')
        save(out/'selection.json', {'included': sorted(set(args.dataset_ids)), 'excluded': sorted(set(available)-set(args.dataset_ids))})
        records = [available[key] for key in sorted(set(args.dataset_ids))]
    if args.run_dials:
        save(out/'failures.json', failures)
    if args.compare:
        tools.run('ccp4-python', [source_dir/'compare_crystals.py', out, args.ccp4, args.reso_cchalf, args.cell_length_tolerance, args.cell_angle_tolerance, args.min_pair_cc, args.min_common, *(args.dataset_ids or [])], out)
    if args.run_mda:
        assemble(tools, records, out, args)
    elif not records:
        raise RuntimeError('All crystals failed; see failures.json')


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError) as error:
        sys.exit(str(error))
