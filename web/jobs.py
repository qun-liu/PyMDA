"""Validated job specifications, isolated outputs and a durable SQLite queue."""
from dataclasses import dataclass
from contextlib import contextmanager
from datetime import datetime, timezone
import json
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
        raise ValueError('Directory is outside the configured data or output roots')
    return path


class JobSpec(BaseModel):
    model_config = ConfigDict(extra='forbid')
    mode: str = 'process'
    input_dir: str = ''
    master_files: list[str] = Field(default_factory=list)
    source_dir: str | None = None
    dataset_ids: list[str] = Field(default_factory=list)
    thread: int = Field(default=2, ge=1, le=2)
    merge_workers: int = Field(default=2, ge=1, le=4)
    wedges: int = Field(default=4, ge=1, le=20)
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
        if value not in {'process', 'assemble', 'both'}:
            raise ValueError('Choose process, assemble or both')
        return value

    @field_validator('prefix_policy', 'crystal_identity')
    @classmethod
    def choice_valid(cls, value, info):
        choices = {'prefix_policy': {'all', 'both', 'selected', 'cchalf', 'full'},
                   'crystal_identity': {'unknown', 'repeat', 'independent'}}
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
    return sorted(Path(source).glob('crystals/*/manifest.json'))


def validate_spec(spec, settings):
    spec = spec.model_copy(deep=True)
    if spec.mode in {'process', 'both'}:
        directory = contained(spec.input_dir or settings.input_root, [settings.input_root])
        if not directory.is_dir():
            raise ValueError('Input directory does not exist')
        available = set(discover_masters(directory))
        selected = [contained(p, [settings.input_root]) for p in spec.master_files]
        if not selected or any(p not in available for p in selected):
            raise ValueError('Select existing rotation master files from the input directory')
        spec.master_files = list(map(str, sorted(set(selected))))
        spec.input_dir = str(directory)
        if spec.mode == 'both' and len(selected) < 2:
            raise ValueError('Assembly requires at least two datasets; choose process for a single scan')
        if spec.mode == 'both' and len(selected) > 10:
            raise ValueError('Select at most ten datasets for the exhaustive assembly search')
    else:
        source = contained(spec.source_dir or settings.output_root, [settings.output_root, settings.jobs_root])
        manifests = {p.parent.name: p for p in source_manifests(source)}
        ids = sorted(set(spec.dataset_ids or manifests))
        if len(ids) < 2 or len(ids) > 10 or not set(ids).issubset(manifests):
            raise ValueError('Select between two and ten processed datasets from the source directory')
        floor = 0.0
        for identifier in ids:
            record = json.loads(manifests[identifier].read_text())
            floor = max(floor, float(record['config']['resolution']))
            for prefix in record['prefixes']:
                path = contained(prefix['mtz'], [source])
                if not path.is_file() or path.stat().st_size == 0:
                    raise ValueError('A source MTZ is missing: %s' % path.name)
        spec.reso = floor
        spec.source_dir = str(source)
        spec.dataset_ids = ids
        spec.input_dir = str(settings.input_root)
    if spec.mode != 'process' and any(r < spec.reso for r in spec.resolutions):
        raise ValueError('Merging cutoffs cannot exceed the integrated resolution (smaller Angstrom values)')
    if spec.reso_cchalf < spec.reso:
        raise ValueError('Prefix assessment cutoff must be at least the integration cutoff')
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
        record['output_dir'] = str(self.settings.jobs_root / record['id'] / 'output')
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
    for identifier in spec.dataset_ids:
        original = source / 'crystals' / identifier / 'manifest.json'
        record = json.loads(original.read_text())
        mappings = {}
        for prefix in record['prefixes']:
            old = contained(prefix['mtz'], [source])
            target = output / old.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_suffix('.copying')
            shutil.copy2(old, temporary)
            temporary.replace(target)
            mappings[str(old)] = str(target)
            prefix['mtz'] = str(target)
        record['selected']['mtz'] = mappings[str(Path(record['selected']['mtz']).resolve())]
        target_manifest = output / 'crystals' / identifier / 'manifest.json'
        temporary = target_manifest.with_suffix('.tmp')
        temporary.write_text(json.dumps(record, indent=2))
        temporary.replace(target_manifest)


def command_for(spec, output, settings):
    command = [settings.python, '-u', str(settings.pipeline), '--ccp4', str(settings.ccp4),
               '--hdf5_path', spec.input_dir, '--dir_output', str(output),
               '--thread', str(spec.thread), '--merge_workers', str(spec.merge_workers),
               '--wedges', str(spec.wedges), '--reso', str(spec.reso),
               '--reso_cchalf', str(spec.reso_cchalf), '--prefix_policy', spec.prefix_policy,
               '--min_completeness', str(spec.min_completeness),
               '--min_outer_completeness', str(spec.min_outer_completeness),
               '--min_cc_half', str(spec.min_cc_half), '--min_isigi', str(spec.min_isigi),
               '--crystal_identity', spec.crystal_identity, '--resolutions', *map(str, spec.resolutions)]
    if spec.mode in {'process', 'both'}:
        command += ['--run_dials', '--master_files', *spec.master_files]
    if spec.mode in {'assemble', 'both'}:
        command += ['--run_mda']
    for field in ('frames', 'spg', 'unit_cell', 'nres'):
        value = getattr(spec, field)
        if value is not None:
            command += ['--' + field, str(value)]
    if spec.fix_detector:
        command.append('--fix_detector')
    return command
