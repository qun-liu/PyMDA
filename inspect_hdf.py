"""Read lightweight Eiger metadata with the CCP4 Python/HDF5 runtime."""
import json
from pathlib import Path
import sys
import h5py

records = []
for name in sys.argv[1:]:
    path = Path(name)
    record = {'path': str(path), 'name': path.name, 'ready': False, 'frames': None, 'companions': [], 'error': ''}
    try:
        with h5py.File(path, 'r') as handle:
            if 'entry/sample/goniometer/omega' in handle:
                record['frames'] = int(handle['entry/sample/goniometer/omega'].size)
            data = handle['entry/data']
            for key in data:
                link = data.get(key, getlink=True)
                if isinstance(link, h5py.ExternalLink):
                    companion = Path(link.filename)
                    if not companion.is_absolute():
                        companion = path.parent / companion
                    if not companion.is_file() and (path.parent / companion.name).is_file():
                        companion = path.parent / companion.name
                    record['companions'].append({'name': companion.name, 'exists': companion.is_file()})
            record['ready'] = all(x['exists'] for x in record['companions'])
            if not record['ready']:
                record['error'] = 'Missing companion diffraction data'
    except Exception as error:
        record['error'] = str(error)
    records.append(record)
print(json.dumps(records))
