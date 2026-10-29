"""Explicit real-data smoke test; requires the localhost server and CCP4 tools."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import httpx


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--url', default='http://127.0.0.1:8000')
    args = parser.parse_args()
    source = Path.home()/'.codex/pymda_output'
    manifests = sorted(source.glob('crystals/*/manifest.json'))
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in manifests}
    with httpx.Client(base_url=args.url, timeout=60) as client:
        def get(path):
            r = client.get(path); r.raise_for_status(); return r
        def submit(spec):
            r = client.post('/api/jobs', json=spec); r.raise_for_status()
            job = r.json(); print('Submitted', job['id'], spec['mode'], flush=True); return job['id']
        data = get('/api/datasets').json()
        assert len(data) == 4 and all(d['ready'] for d in data)
        html = get('/').text
        assert '<script' not in html.lower()
        assembly = submit(dict(mode='assemble', source_dir=str(source), dataset_ids=[p.parent.name for p in manifests[:2]], resolutions=[2.5], prefix_policy='full', merge_workers=1))
        process = submit(dict(mode='process', master_files=[data[0]['path']], frames=180, wedges=2, thread=1))
        pending = {assembly, process}
        deadline = time.monotonic()+1800
        previous = {}
        while pending and time.monotonic() < deadline:
            for identifier in list(pending):
                job = get('/api/jobs/'+identifier).json()
                if previous.get(identifier) != job['status']:
                    print(identifier, job['status'], job.get('error'), flush=True)
                    previous[identifier] = job['status']
                if job['status'] in {'succeeded','partial','failed','cancelled','interrupted'}:
                    pending.remove(identifier)
                    print(json.dumps(job, indent=2), flush=True)
                    assert job['status'] == 'succeeded', job
                    page = get('/jobs/'+identifier+'?refresh=0').text
                    assert '<script' not in page.lower()
                    if identifier == assembly:
                        mtz = get('/jobs/'+identifier+'/file?path=best/final.mtz').content
                        assert mtz[:4] == b'MTZ '
                        assert job['statistics']
            if pending:
                time.sleep(5)
        assert not pending, 'Timed out; jobs remain visible in the interface'
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p,h in hashes.items())
    print('Real-data web smoke tests passed; source manifests unchanged.', flush=True)

if __name__ == '__main__':
    main()
