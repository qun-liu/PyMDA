"""Explicit real-data smoke test; requires the localhost server and CCP4 tools."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import httpx
import re
from html import unescape


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--url', default='http://127.0.0.1:8000')
    parser.add_argument('--source', default=str(Path.home()/'.codex/pymda_output'))
    parser.add_argument('--input-dir', default=str(Path.home()/'.codex/pymda_data'))
    parser.add_argument('--skip-process', action='store_true', help='Reuse the separately validated raw-image processing')
    parser.add_argument('--output-dir', help='Verify jobs write to this selected existing output folder')
    parser.add_argument('--resolution',type=float,default=2.5,help='Resolution for assembly and comparison (Angstrom)')
    args = parser.parse_args()
    source = Path(args.source).resolve()
    with httpx.Client(base_url=args.url, timeout=60) as client:
        def get(path):
            r = client.get(path); r.raise_for_status(); return r
        def submit(spec):
            if args.output_dir:
                spec['output_dir'] = str(Path(args.output_dir).expanduser().resolve())
            r = client.post('/api/jobs', json=spec); r.raise_for_status()
            job = r.json(); print('Submitted', job['id'], spec['mode'], flush=True); return job['id']
        catalog=client.get('/api/processed-datasets',params={'source_dir':str(source)})
        catalog.raise_for_status(); catalog=catalog.json()
        selected=[entry for entry in catalog if entry['ready'] and entry['default_selected']][:2]
        assert selected, 'This smoke test requires completed scans in the selected source folder'
        assert args.resolution>=max(entry['integration_resolution'] for entry in selected), 'Choose a cutoff supported by the selected integrations'
        manifests={entry['id']:Path(entry['source_dir'])/'crystals'/entry['dataset_id']/'manifest.json' for entry in catalog}
        hashes={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in manifests.values()}
        for entry in selected:
            record=json.loads(manifests[entry['id']].read_text())
            for export in record['prefixes'] + [w for w in record.get('scan_wedges',[]) if w.get('mtz')]:
                path=Path(export['mtz'])
                hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        html=client.get('/',params={'source_dir':str(source)}); html.raise_for_status(); html=html.text
        assert 'Find processed datasets in this output folder' in html
        assert 'Open this result directory' in html
        checkboxes=re.findall(r'<input[^>]*name="dataset_ids"[^>]*>',html)
        assert sum('checked' in tag for tag in checkboxes)==sum(entry['default_selected'] for entry in catalog)
        for entry in selected:
            assert entry['id'] in html and entry['source_dir'] in html
        data = client.get('/api/datasets',params={'input_dir':args.input_dir}); data.raise_for_status(); data=data.json()
        assert data and all(d['ready'] for d in data)
        html = get('/').text
        assert '<script' not in html.lower()
        assert 'Consecutive wedges per scan' in html
        assert 'Dataset subset search' in html
        if args.output_dir:
            directory=Path(args.output_dir).expanduser().resolve()
            response=client.get('/browse',params={'kind':'output_dir','path':str(directory),'input_dir':args.input_dir,'source_dir':str(source)})
            response.raise_for_status()
            chosen=unescape(re.search(r'href="([^"]+)">Use this directory</a>',response.text).group(1))
            response=client.get(chosen); response.raise_for_status()
            assert str(directory) in get('/').text
        identifiers=[entry['id'] for entry in selected]
        print('Discovered',len(catalog),'processing results; selected',identifiers,flush=True)
        comparison = submit(dict(mode='compare', source_dir=str(source), dataset_ids=identifiers, resolutions=[args.resolution],reso_cchalf=args.resolution))
        single = submit(dict(mode='assemble', source_dir=str(source), dataset_ids=identifiers[:1], resolutions=[args.resolution],reso_cchalf=args.resolution, prefix_policy='full', merge_workers=1))
        expected={comparison:identifiers,single:identifiers[:1]}
        assembly=None
        if len(selected)>1:
            assembly=submit(dict(mode='assemble', source_dir=str(source), dataset_ids=identifiers, resolutions=[args.resolution],reso_cchalf=args.resolution,prefix_policy='full', merge_workers=1))
            expected[assembly]=identifiers
        pending = set(expected)
        if not args.skip_process:
            pending.add(submit(dict(mode='process', input_dir=args.input_dir, master_files=[data[0]['path']], frames=180, wedges=2, thread=1)))
        deadline = time.monotonic()+1800
        previous = {}
        while pending and time.monotonic() < deadline:
            for identifier in list(pending):
                job = get('/api/jobs/'+identifier).json()
                if args.output_dir:
                    assert Path(job['output_dir']).is_relative_to(directory)
                if previous.get(identifier) != job['status']:
                    print(identifier, job['status'], job.get('error'), flush=True)
                    previous[identifier] = job['status']
                if job['status'] in {'succeeded','partial','failed','cancelled','interrupted'}:
                    pending.remove(identifier)
                    print(json.dumps(job, indent=2), flush=True)
                    if identifier == single and job['status'] == 'failed':
                        trials = get('/jobs/'+identifier+'/file?path=trials.json').json()
                        assert trials and all(not t['eligible'] for t in trials)
                        assert client.get('/jobs/'+identifier+'/file?path=best/final.mtz').status_code == 404
                        continue
                    assert job['status'] == 'succeeded', job
                    page = get('/jobs/'+identifier+'?refresh=0').text
                    assert '<script' not in page.lower()
                    if identifier in expected:
                        selection=get('/jobs/'+identifier+'/file?path=selection.json').json()
                        assert set(selection['included'])==set(expected[identifier])
                        assert selection['source_manifests']=={key:str(manifests[key]) for key in expected[identifier]}
                        for key in expected[identifier]:
                            record=get('/jobs/'+identifier+'/file?path=crystals/'+key+'/manifest.json').json()
                            assert record['dataset_id']==key and record['source_manifest']==str(manifests[key])
                            assert Path(record['selected']['mtz']).is_relative_to(Path(job['output_dir'])/'crystals'/key)
                    if identifier == comparison:
                        report = get('/jobs/'+identifier+'/file?path=comparison.json').json()
                        assert len(report['datasets']) == len(selected) and len(report['pairs']) == len(selected)*(len(selected)-1)//2
                        assert {entry['id'] for entry in report['datasets']}==set(identifiers)
                        assert all(pair['correlation'] is not None for pair in report['pairs'])
                        assert 'Select group' in page
                        group_link=unescape(re.search(r'href="([^"]+)">Select group 1',page).group(1))
                        group=get(group_link).text
                        group_choices=re.findall(r'<input[^>]*name="dataset_ids"[^>]*>',group)
                        checked={unescape(re.search(r'value="([^"]+)"',tag).group(1)) for tag in group_choices if 'checked' in tag}
                        assert checked==set(report['groups'][0])
                    if identifier in {assembly,single}:
                        mtz = get('/jobs/'+identifier+'/file?path=best/final.mtz').content
                        assert mtz[:4] == b'MTZ '
                        assert job['statistics']
                        plan = get('/jobs/'+identifier+'/file?path=search.json').json()
                        assert plan['dataset_count'] in {1,2}
                        assert plan['candidate_count'] >= 1
                        assert set(plan['dataset_subsets'][0])==set(expected[identifier])
                        record = get('/jobs/'+identifier+'/file?path=crystals/'+identifiers[0]+'/manifest.json').json()
                        for wedge in record.get('scan_wedges', []):
                            assert f"{wedge['first_frame']}–{wedge['last_frame']}" in page
                            relative = Path(wedge['mtz']).relative_to(Path(job['output_dir']))
                            result = client.get('/jobs/'+identifier+'/file',params={'path':str(relative)})
                            result.raise_for_status()
                            assert result.content[:4] == b'MTZ '
            if pending:
                time.sleep(5)
        assert not pending, 'Timed out; jobs remain visible in the interface'
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p,h in hashes.items())
    print('Real-data web smoke tests passed; source manifests and exported MTZs unchanged.', flush=True)

if __name__ == '__main__':
    main()
