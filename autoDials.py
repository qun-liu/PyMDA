"""Checked subprocess execution for current CCP4/DIALS installations."""
import json
import os
from pathlib import Path
import shlex
import subprocess


class Toolchain:
    def __init__(self, root):
        root = Path(root).expanduser().resolve()
        if not (root / 'bin/ccp4.setup-sh').is_file():
            candidates = sorted(root.glob('ccp4-*/bin/ccp4.setup-sh'))
            if not candidates:
                raise ValueError('No CCP4 setup script below %s' % root)
            root = candidates[-1].parent.parent
        self.root = root
        setup = root / 'bin/ccp4.setup-sh'
        result = subprocess.run(['bash', '-c', 'source "$1" >/dev/null && env -0', 'pymda', str(setup)], check=True, stdout=subprocess.PIPE)
        self.env = dict(item.split('=', 1) for item in result.stdout.decode().split('\0') if '=' in item)
        self.env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg')
        for name in ('dials.import', 'dials.find_spots', 'dials.index', 'dials.refine', 'dials.integrate', 'dials.slice_sequence', 'dials.export', 'pointless', 'aimless', 'ctruncate', 'freerflag', 'ccp4-python'):
            if not os.access(root / 'bin' / name, os.X_OK):
                raise ValueError('Missing executable: %s' % (root / 'bin' / name))

    def run(self, name, args, directory, stdin=None, expected=()):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        command = [str(self.root / 'bin' / name)] + [str(x) for x in args]
        log = directory / (name + '.console.log')
        input_files = {}
        file_arguments = list(args)
        if stdin:
            for line in stdin.splitlines():
                fields = shlex.split(line)
                if len(fields) == 2 and fields[0].upper() == 'HKLIN':
                    file_arguments.append(fields[1])
        for arg in file_arguments:
            candidate = Path(str(arg))
            if not candidate.is_absolute():
                candidate = directory / candidate
            if candidate.is_file() and candidate.name not in expected:
                stat = candidate.stat()
                input_files[str(candidate.resolve())] = [stat.st_size, stat.st_mtime_ns]
        specification = {'argv': command, 'stdin': stdin, 'input_files': input_files}
        command_file = directory / (name + '.command.json')
        success_file = directory / (name + '.success.json')
        if expected and success_file.is_file() and json.loads(success_file.read_text()) == specification and all((directory / f).is_file() and (directory / f).stat().st_size for f in expected):
            return log
        if success_file.exists():
            success_file.unlink()
        command_file.write_text(json.dumps(specification, indent=2))
        for filename in expected:
            path = directory / filename
            if path.exists():
                path.unlink()
        print('%s: %s' % (directory.name, shlex.join(command)), flush=True)
        with log.open('w') as out:
            result = subprocess.run(command, cwd=directory, env=self.env, input=stdin, text=True, stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError('%s exited %s; see %s' % (name, result.returncode, log))
        for filename in expected:
            path = directory / filename
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError('%s did not produce %s; see %s' % (name, path, log))
        success_file.write_text(json.dumps(specification, indent=2))
        return log
