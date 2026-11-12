"""Validated job specifications, isolated outputs and a durable SQLite queue."""
from dataclasses import dataclass
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import uuid

from pydantic import BaseModel, ConfigDict, Field, field_validator
from helpers import discover_masters

ROOT = Path(__file__).resolve().parent.parent
TERMINAL = {'succeeded', 'partial', 'failed', 'cancelled', 'interrupted'}
RESUMABLE = {'partial', 'failed', 'cancelled', 'interrupted'}


def now():
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Settings:
    input_root: Path
    output_root: Path
    jobs_root: Path
    ccp4: Path = Path('/home/aa1/ccp47')
    pipeline: Path = ROOT / 'pymda.py'
    python: str = sys.executable
    testing: bool = False

    @classmethod
    def from_env(cls):
        output = Path(os.getenv('PYMDA_OUTPUT_ROOT', '~/.codex/pymda_output')).expanduser().resolve()
        return cls(Path(os.getenv('PYMDA_INPUT_ROOT', '~/.codex/pymda_data')).expanduser().resolve(),
                   output, Path(os.getenv('PYMDA_JOBS_ROOT', str(output / 'jobs'))).expanduser().resolve(),
                   Path(os.getenv('PYMDA_CCP4', '/home/aa1/ccp47')).expanduser().resolve())

    @property
    def database(self):
        return self.jobs_root / 'queue.sqlite3'


def contained(path, roots):
    path = Path(path).expanduser().resolve()
    if not any(path == Path(r).resolve() or Path(r).resolve() in path.parents for r in roots):
        raise ValueError('File is outside the result directory')
    return path


def local_directory(path, writable=False):
    directory = Path(path).expanduser().resolve()
    if not directory.is_dir():
        if directory.exists():
            raise ValueError('Path is not a directory: %s' % directory)
        raise ValueError('Directory does not exist: %s' % directory)
    if not os.access(directory, (os.W_OK if writable else os.R_OK) | os.X_OK):
        message = 'Output directory is not writable' if writable else 'Directory is not readable'
        raise ValueError('%s: %s' % (message, directory))
    return directory


def default_directory(*paths, writable=False):
    """Use an accessible starting folder when saved/default locations disappear."""
    for path in (*paths, Path.home(), Path('/')):
        if not path:
            continue
        try:
            return local_directory(path, writable=writable)
        except (ValueError, OSError):
            pass
    raise ValueError('No accessible %s directory found' % ('output' if writable else 'starting'))


def output_directory(path, create=False):
    """Validate an output folder or its existing parent; create only on request."""
    directory = Path(path).expanduser().resolve()
    parent = directory
    while not parent.exists():
        parent = parent.parent
    local_directory(parent, writable=True)
    if create:
        directory.mkdir(parents=True, exist_ok=True)
        local_directory(directory, writable=True)
    return directory


class JobSpec(BaseModel):
    model_config = ConfigDict(extra='forbid')
    cell_length_tolerance: float = Field(default=2, gt=0, allow_inf_nan=False)
    cell_angle_tolerance: float = Field(default=2, gt=0, allow_inf_nan=False)
    min_pair_cc: float = Field(default=0.9, ge=-1, le=1, allow_inf_nan=False)
    min_common: int = Field(default=100, ge=3)
    output_dir: str = ''
    mode: str = 'process'
    input_dir: str = ''
    master_files: list[str] = Field(default_factory=list)
    source_dir: str | None = None
    dataset_ids: list[str] = Field(default_factory=list)
    thread: int = Field(default=2, ge=1, le=2)
    merge_workers: int = Field(default=2, ge=1, le=4)
    wedges: int = Field(default=4, ge=1, le=20)
    split_wedges: int = Field(default=1, ge=1, le=100)
    subset_search: str = 'auto'
    max_subsets: int = Field(default=64, ge=1, le=4096)
    frames: int | None = Field(default=None, ge=1)
    reso: float = Field(default=1.8, gt=0, allow_inf_nan=False)
    reso_cchalf: float = Field(default=2.5, gt=0, allow_inf_nan=False)
    resolutions: list[float] = Field(default_factory=lambda: [2.0, 2.2, 2.26, 2.3, 2.5, 3.0], min_length=1, max_length=30)
    spg: str | None = None
    unit_cell: str | None = None
    fix_detector: bool = False
    min_completeness: float = Field(default=90, ge=0, le=100, allow_inf_nan=False)
    min_outer_completeness: float = Field(default=80, ge=0, le=100, allow_inf_nan=False)
    min_cc_half: float = Field(default=0.3, ge=0, le=1, allow_inf_nan=False)
    min_isigi: float = Field(default=1, ge=0, allow_inf_nan=False)
    prefix_policy: str = 'all'
    nres: int | None = Field(default=None, ge=1)
    crystal_identity: str = 'unknown'

    @field_validator('mode')
    @classmethod
    def mode_valid(cls, value):
        if value not in {'process', 'assemble', 'both', 'compare'}:
            raise ValueError('Choose process, compare, assemble or both')
        return value

    @field_validator('prefix_policy', 'crystal_identity', 'subset_search')
    @classmethod
    def choice_valid(cls, value, info):
        choices = {'prefix_policy': {'all', 'both', 'selected', 'cchalf', 'full'},
                   'crystal_identity': {'unknown', 'repeat', 'independent'},
                   'subset_search': {'auto', 'exhaustive', 'bounded'}}
        if value not in choices[info.field_name]:
            raise ValueError('Invalid selection')
        return value

    @field_validator('resolutions')
    @classmethod
    def resolutions_valid(cls, values):
        if any(not math.isfinite(v) or v <= 0 or abs(v - round(v, 2)) > 1e-8 for v in values):
            raise ValueError('Cutoffs must be positive finite values with at most two decimal places')
        return sorted(set(values))

    @field_validator('spg')
    @classmethod
    def spg_valid(cls, value):
        if value and not re.fullmatch(r'[A-Za-z0-9 /:_-]{1,50}', value):
            raise ValueError('Invalid space-group symbol')
        return value or None

    @field_validator('unit_cell')
    @classmethod
    def cell_valid(cls, value):
        if not value:
            return None
        try:
            numbers = [float(x) for x in value.replace(',', ' ').split()]
        except ValueError:
            raise ValueError('Enter six numeric unit-cell parameters') from None
        if len(numbers) != 6 or not all(math.isfinite(x) for x in numbers):
            raise ValueError('Enter six finite unit-cell parameters')
        if min(numbers[:3]) <= 0 or not all(0 < x < 180 for x in numbers[3:]):
            raise ValueError('Lengths must be positive and angles between 0 and 180 degrees')
        a, b, c = [math.cos(math.radians(x)) for x in numbers[3:]]
        if 1 - a*a - b*b - c*c + 2*a*b*c <= 0:
            raise ValueError('Unit-cell angles must define a positive volume')
        return ','.join(map(str, numbers))


def source_manifests(source):
    """Find direct results and results nested under job/output directories."""
    root = Path(source).resolve()
    if root.name == 'crystals':
        return sorted(root.glob('*/manifest.json'))
    for directory in (root, *root.parents):
        if directory.parent.name == 'crystals':
            manifest = directory/'manifest.json'
            return [manifest] if manifest.is_file() else []
    return sorted(root.glob('**/crystals/*/manifest.json'))


def source_manifest_map(source):
    root = Path(source).resolve()
    manifests = {}
    for path in source_manifests(root):
        identifier = path.parent.name
        result_root = path.parent.parent.parent
        if result_root != root and result_root not in root.parents:
            relative = path.relative_to(root)
            identifier = identifier[:120] + '__' + hashlib.sha256(relative.as_posix().encode()).hexdigest()[:16]
        manifests[identifier] = path
    return manifests


def source_mtz(manifest, reference):
    """Resolve an export inside its crystal directory, including moved results."""
    crystal = Path(manifest).parent
    requested = Path(reference).expanduser()
    candidates = [requested] if requested.is_absolute() else [crystal/requested, crystal.parent.parent/requested]
    for index, part in enumerate(requested.parts[:-1]):
        if part == 'crystals' and requested.parts[index+1] == crystal.name:
            candidates.append(crystal.joinpath(*requested.parts[index+2:]))
    for candidate in candidates:
        try:
            path = contained(candidate, [crystal])
            if path.is_file() and path.stat().st_size and os.access(path, os.R_OK):
                return path
        except (ValueError,OSError):
            pass
    raise ValueError('Missing exported MTZ in %s: %s' % (crystal, requested.name))


def source_catalog(source):
    entries, preferred = [], {}
    for identifier, path in source_manifest_map(source).items():
        entry = {'id': identifier, 'dataset_id': path.parent.name, 'manifest': path,
                 'source_dir': str(path.parent.parent.parent), 'name': path.parent.name,
                 'frame': None, 'integration_resolution': None, 'ready': False,
                 'default_selected': False, 'error': '', 'record': None}
        try:
            record = json.loads(path.read_text())
            entry['name'] = Path(record['config']['master']).name
            resolution = float(record['config']['resolution'])
            if not math.isfinite(resolution) or resolution <= 0 or not record['prefixes']:
                raise ValueError('Incomplete processing manifest')
            entry['integration_resolution'] = resolution
            exports = set()
            for prefix in record['prefixes'] + [record['selected']]:
                if not isinstance(prefix['last_frame'], int) or prefix['last_frame'] <= 0:
                    raise ValueError('Invalid last frame in processing result')
                for key in ('resolution','isigi','cc_half','completeness'):
                    values = prefix['statistics'][key]
                    if not isinstance(values, list) or len(values) < 3 or not all(math.isfinite(float(v)) for v in values):
                        raise ValueError('Incomplete prefix statistics')
            entry['frame'] = record['selected']['last_frame']
            for prefix in record['prefixes']:
                exports.add(source_mtz(path, prefix['mtz']))
            if source_mtz(path, record['selected']['mtz']) not in exports:
                raise ValueError('Selected MTZ is not one of the processed exports')
            for wedge in record.get('scan_wedges', []):
                if wedge.get('mtz'):
                    source_mtz(path, wedge['mtz'])
            score = (record['config'].get('frames') is None,
                     max(p['last_frame'] for p in record['prefixes']),
                     -resolution, path.stat().st_mtime_ns)
            entry.update(ready=True, record=record)
            master = record['config']['master']
            if master not in preferred or score > preferred[master][0]:
                preferred[master] = (score, entry)
        except (ValueError,OSError,KeyError,TypeError) as error:
            entry['error'] = 'Incomplete or unavailable processing result: %s' % error
        entries.append(entry)
    for score, entry in preferred.values():
        entry['default_selected'] = True
    return entries


def validate_spec(spec, settings):
    spec = spec.model_copy(deep=True)
    if spec.mode in {'process', 'both'}:
        if spec.frames is not None and spec.split_wedges > spec.frames:
            raise ValueError('Consecutive wedge count cannot exceed the pilot frame limit')
        directory = local_directory(spec.input_dir) if spec.input_dir else default_directory(
            os.getenv('PYMDA_START_INPUT_DIR'), settings.input_root, ROOT)
        available = set(discover_masters(directory))
        selected = [Path(p).expanduser().resolve() for p in spec.master_files]
        if not selected or any(p not in available for p in selected):
            raise ValueError('Select existing rotation master files from the input directory')
        for path in selected:
            if not os.access(path, os.R_OK):
                raise ValueError('Master file is not readable: %s' % path)
        spec.master_files = list(map(str, sorted(set(selected))))
        spec.input_dir = str(directory)
        if spec.mode == 'both' and spec.subset_search == 'exhaustive' and len(spec.master_files) > 10:
            raise ValueError('Use auto or bounded search for more than ten datasets')
    else:
        source = local_directory(spec.source_dir) if spec.source_dir else default_directory(
            os.getenv('PYMDA_START_SOURCE_DIR'), settings.output_root, settings.jobs_root, ROOT)
        entries = {entry['id']: entry for entry in source_catalog(source)}
        ids = sorted(set(spec.dataset_ids or [key for key,entry in entries.items() if entry['default_selected']]))
        if not ids or not set(ids).issubset(entries):
            raise ValueError('Select one or more processed datasets from the source directory')
        if spec.mode == 'assemble' and spec.subset_search == 'exhaustive' and len(ids) > 10:
            raise ValueError('Use auto or bounded search for more than ten datasets')
        floor = 0.0
        for identifier in ids:
            entry = entries[identifier]
            if not entry['ready']:
                raise ValueError(entry['error'])
            record = entry['record']
            floor = max(floor, float(record['config']['resolution']))
        spec.reso = floor
        spec.source_dir = str(source)
        spec.dataset_ids = ids
        spec.input_dir = str(settings.input_root)
    if spec.mode in {'assemble','both'} and any(r < spec.reso for r in spec.resolutions):
        raise ValueError('Merging cutoffs cannot exceed the integrated resolution (smaller Angstrom values)')
    if spec.reso_cchalf < spec.reso:
        raise ValueError('Prefix assessment cutoff must be at least the integration cutoff')
    if spec.output_dir:
        spec.output_dir = str(output_directory(spec.output_dir))
    return spec


class JobStore:
    def __init__(self, settings):
        self.settings = settings

    @contextmanager
    def connect(self):
        connection = sqlite3.connect(self.settings.database, timeout=15)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def initialize(self):
        self.settings.jobs_root.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.execute('PRAGMA journal_mode=WAL')
            db.execute('''CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY, status TEXT NOT NULL, params TEXT NOT NULL,
                created TEXT NOT NULL, started TEXT, finished TEXT,
                pid INTEGER, cancel_requested INTEGER NOT NULL DEFAULT 0,
                message TEXT NOT NULL DEFAULT '', stage TEXT NOT NULL DEFAULT 'Queued')''')

    def unpack(self, row):
        if row is None:
            return None
        record = dict(row)
        record['params'] = json.loads(record['params'])
        record['output_dir'] = str(Path(record['params']['output_dir']) / record['id'] / 'output') if record['params'].get('output_dir') else str(self.settings.jobs_root / record['id'] / 'output')
        return record

    def get(self, identifier):
        with self.connect() as db:
            return self.unpack(db.execute('SELECT * FROM jobs WHERE id=?', [identifier]).fetchone())

    def list(self, limit=50):
        with self.connect() as db:
            return [self.unpack(row) for row in db.execute('SELECT * FROM jobs ORDER BY created DESC LIMIT ?', [limit])]

    def create(self, spec):
        identifier = uuid.uuid4().hex
        with self.connect() as db:
            db.execute('INSERT INTO jobs(id,status,params,created) VALUES(?,?,?,?)',
                       [identifier, 'queued', spec.model_dump_json(), now()])
        return self.get(identifier)

    def claim(self):
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            if db.execute("SELECT 1 FROM jobs WHERE status='running' LIMIT 1").fetchone():
                return None
            row = db.execute("SELECT * FROM jobs WHERE status='queued' ORDER BY created LIMIT 1").fetchone()
            if row is None:
                return None
            db.execute("UPDATE jobs SET status='running',started=?,finished=NULL,stage='Preparing inputs' WHERE id=?", [now(), row['id']])
        return self.get(row['id'])

    def update(self, identifier, **fields):
        allowed = {'status', 'started', 'finished', 'pid', 'cancel_requested', 'message', 'stage'}
        if not fields or not set(fields).issubset(allowed):
            raise ValueError('Invalid job update')
        with self.connect() as db:
            db.execute('UPDATE jobs SET ' + ','.join(k + '=?' for k in fields) + ' WHERE id=?', [*fields.values(), identifier])

    def cancel(self, identifier):
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            row = db.execute('SELECT status FROM jobs WHERE id=?', [identifier]).fetchone()
            if row is None:
                raise KeyError(identifier)
            if row['status'] == 'queued':
                db.execute("UPDATE jobs SET status='cancelled',cancel_requested=1,finished=?,stage='Cancelled' WHERE id=?", [now(), identifier])
            elif row['status'] == 'running':
                db.execute("UPDATE jobs SET cancel_requested=1,stage='Cancelling' WHERE id=?", [identifier])
            else:
                raise ValueError('Only queued or running jobs can be cancelled')

    def resume(self, identifier):
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            row = db.execute('SELECT status FROM jobs WHERE id=?', [identifier]).fetchone()
            if row is None:
                raise KeyError(identifier)
            if row['status'] not in RESUMABLE:
                raise ValueError('Only failed, partial, cancelled or interrupted jobs can be resumed')
            db.execute("UPDATE jobs SET status='queued',pid=NULL,cancel_requested=0,finished=NULL,message='',stage='Queued for recovery' WHERE id=?", [identifier])


def prepare_assembly(job, spec):
    """Copy only exported data and remap every MTZ reference into this job."""
    source = Path(spec.source_dir)
    output = Path(job['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    comparison = source/'comparison.json'
    if comparison.is_file():
        shutil.copy2(comparison, output/'comparison.json')
    manifests = source_manifest_map(source)
    (output/'selection.json').write_text(json.dumps({'included':spec.dataset_ids,
        'excluded':[key for key in manifests if key not in spec.dataset_ids], 'source':str(source),
        'source_manifests':{key:str(manifests[key]) for key in spec.dataset_ids}},indent=2))
    for identifier in spec.dataset_ids:
        original = manifests[identifier]
        record = json.loads(original.read_text())
        selected = source_mtz(original, record['selected']['mtz'])
        mappings = {}
        for prefix in record['prefixes'] + [w for w in record.get('scan_wedges', []) if w.get('mtz')]:
            old = source_mtz(original, prefix['mtz'])
            target = output/'crystals'/identifier/old.relative_to(original.parent.resolve())
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_suffix('.copying')
            shutil.copy2(old, temporary)
            temporary.replace(target)
            mappings[str(old)] = str(target)
            prefix['mtz'] = str(target)
        record['selected']['mtz'] = mappings[str(selected)]
        record['dataset_id'] = identifier
        record['source_manifest'] = str(original)
        target_manifest = output / 'crystals' / identifier / 'manifest.json'
        temporary = target_manifest.with_suffix('.tmp')
        temporary.write_text(json.dumps(record, indent=2))
        temporary.replace(target_manifest)


def command_for(spec, output, settings):
    command = [settings.python, '-u', str(settings.pipeline), '--ccp4', str(settings.ccp4),
               '--hdf5_path', spec.input_dir, '--dir_output', str(output),
               '--thread', str(spec.thread), '--merge_workers', str(spec.merge_workers),
               '--wedges', str(spec.wedges), '--reso', str(spec.reso),
               '--subset_search', spec.subset_search, '--max_subsets', str(spec.max_subsets),
               '--reso_cchalf', str(spec.reso_cchalf), '--prefix_policy', spec.prefix_policy,
               '--min_completeness', str(spec.min_completeness),
               '--min_outer_completeness', str(spec.min_outer_completeness),
               '--min_cc_half', str(spec.min_cc_half), '--min_isigi', str(spec.min_isigi),
               '--crystal_identity', spec.crystal_identity, '--resolutions', *map(str, spec.resolutions)]
    if spec.mode in {'process', 'both'}:
        command += ['--split_wedges', str(spec.split_wedges), '--run_dials', '--master_files', *spec.master_files]
    if spec.mode == 'compare':
        command += ['--compare', '--cell_length_tolerance', str(spec.cell_length_tolerance), '--cell_angle_tolerance', str(spec.cell_angle_tolerance), '--min_pair_cc', str(spec.min_pair_cc), '--min_common', str(spec.min_common)]
    if spec.mode in {'assemble', 'both'}:
        command += ['--run_mda']
    for field in ('frames', 'spg', 'unit_cell', 'nres'):
        value = getattr(spec, field)
        if value is not None:
            command += ['--' + field, str(value)]
    if spec.fix_detector:
        command.append('--fix_detector')
    return command
