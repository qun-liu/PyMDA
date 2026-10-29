#!/usr/bin/env python3
"""A single persistent worker. Scientific subprocesses have their own process group."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from web.jobs import JobSpec, JobStore, Settings, now, prepare_assembly, command_for, validate_spec

STOP = False


def stop(signum, frame):
    global STOP
    STOP = True


def owned_groups(identifier):
    """Identify descendants by their inherited UUID marker, avoiding reused PIDs."""
    marker = ('PYMDA_WEB_JOB_ID=' + identifier).encode()
    groups = set()
    for directory in Path('/proc').glob('[0-9]*'):
        try:
            if marker not in (directory / 'environ').read_bytes().split(b'\0'):
                continue
            fields = (directory / 'stat').read_text().rsplit(')', 1)[1].split()
            groups.add(int(fields[2]))
        except (OSError, ValueError, IndexError):
            continue
    return groups


def terminate_job(identifier):
    groups = owned_groups(identifier)
    for group in groups:
        if group == os.getpgrp():
            continue
        try:
            os.killpg(group, signal.SIGTERM)
        except ProcessLookupError:
            pass
    if groups:
        time.sleep(1)
    for group in groups:
        if group == os.getpgrp():
            continue
        try:
            os.killpg(group, signal.SIGKILL)
        except ProcessLookupError:
            pass


def read_json(path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def stage_from_log(path):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size - 12000))
            text = stream.read().decode(errors='replace')
        names = {'dials.import': 'Importing images', 'dials.find_spots': 'Finding diffraction spots',
                 'dials.index': 'Indexing', 'dials.refine': 'Refining geometry',
                 'dials.integrate': 'Integrating reflections', 'dials.slice_sequence': 'Exporting frame prefixes',
                 'dials.export': 'Exporting frame prefixes', 'pointless': 'Checking and combining data',
                 'aimless': 'Scaling and comparing merges', 'ctruncate': 'Converting intensities to amplitudes',
                 'freerflag': 'Assigning Free-R flags', 'ccp4-python': 'Plotting results'}
        for line in reversed(text.splitlines()):
            for tool, label in names.items():
                if '/bin/' + tool + ' ' in line:
                    return label
    except OSError:
        pass
    return 'Running PyMDA'


def run_job(store, job):
    process = None
    output = Path(job['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    log = output / 'job.stdout.log'
    try:
        spec = validate_spec(JobSpec(**job['params']), store.settings)
        if spec.mode == 'assemble':
            prepare_assembly(job, spec)
        if STOP or store.get(job['id'])['cancel_requested']:
            status = 'interrupted' if STOP else 'cancelled'
            store.update(job['id'], status=status, finished=now(), stage=status.title())
            return
        command = command_for(spec, output, store.settings)
        (output / 'web_command.json').write_text(json.dumps(command, indent=2))
        environment = dict(os.environ, PYMDA_WEB_JOB_ID=job['id'])
        with log.open('a') as stream:
            stream.write('\n--- Job attempt started %s ---\n' % now())
            stream.flush()
            process = subprocess.Popen(command, cwd=Path(__file__).resolve().parent,
                                       stdout=stream, stderr=subprocess.STDOUT,
                                       env=environment, start_new_session=True)
            store.update(job['id'], pid=process.pid, stage='Running PyMDA')
            while process.poll() is None:
                current = store.get(job['id'])
                if STOP or current['cancel_requested']:
                    terminate_job(job['id'])
                    process.wait(timeout=10)
                    status = 'interrupted' if STOP else 'cancelled'
                    store.update(job['id'], status=status, finished=now(), pid=None,
                                 stage=status.title(), message='Completed stages are retained for recovery.')
                    return
                store.update(job['id'], stage=stage_from_log(log))
                time.sleep(0.25)
        # Catch children left behind if the pipeline itself exited unexpectedly.
        terminate_job(job['id'])
        current = store.get(job['id'])
        if current['cancel_requested']:
            status, message = 'cancelled', 'Completed stages are retained for recovery.'
        elif process.returncode:
            status, message = 'failed', 'PyMDA exited with code %s. See the log for details.' % process.returncode
        elif spec.mode != 'process' and not (output / 'best/final.mtz').is_file():
            status, message = 'failed', 'The pipeline did not produce a final MTZ.'
        elif spec.mode == 'process' and not any(output.glob('crystals/*/manifest.json')):
            status, message = 'failed', 'The pipeline did not produce processed dataset manifests.'
        else:
            failures = read_json(output / 'failures.json', [])
            status = 'partial' if failures else 'succeeded'
            message = '%s dataset(s) failed; inspect the log and resume to retry.' % len(failures) if failures else 'Completed successfully.'
        store.update(job['id'], status=status, message=message, finished=now(), pid=None, stage=status.title())
    except Exception as error:
        if process is not None:
            terminate_job(job['id'])
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
        with log.open('a') as stream:
            stream.write('\nWorker: %s\n' % error)
        store.update(job['id'], status='interrupted' if STOP else 'failed', finished=now(), pid=None,
                     stage='Interrupted' if STOP else 'Failed', message=str(error))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    settings = Settings(**{k: Path(v) if k in {'input_root', 'output_root', 'jobs_root', 'ccp4', 'pipeline'} else v for k, v in config.items()})
    store = JobStore(settings)
    store.initialize()
    with (settings.jobs_root / 'worker.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, stop)
        for job in store.list(100000):
            if job['status'] == 'running':
                terminate_job(job['id'])
                store.update(job['id'], status='interrupted', pid=None, finished=now(), stage='Interrupted',
                             message='Worker restarted. Resume to reuse completed stages.')
        while not STOP:
            job = store.claim()
            if job:
                run_job(store, job)
            else:
                time.sleep(0.25)


if __name__ == '__main__':
    main()
