# Local PyMDA interface

Start from the repository directory:

```bash
source .pymda/bin/activate
python run_web.py
```

Open http://127.0.0.1:8000. Stop the server with Ctrl-C. The interface uses
FastAPI and Python-generated HTML/CSS. Forms and five-second page refreshes
work without JavaScript or Jinja2. There is no authentication; the server binds
to loopback and accepts local requests only.

The workspace recursively discovers diffraction master HDF5 files and checks
frame counts and companion files with CCP4 Python. Select scans to process,
process and assemble, or select existing processed datasets for assembly alone.
Advanced parameters include pilot frame limits, geometry, resolution cutoffs,
quality thresholds and crystal provenance. Use a pilot before a full scan run.

A separate worker runs one queued job at a time. SQLite preserves job state;
each job has an isolated output directory. Job pages show status, selectable
logs, merge statistics, resolution tradeoff plots and downloadable artifacts.
Cancel terminates the job's subprocess group. If the server or worker stops,
unfinished jobs become interrupted. Resume retries the same parameters and
reuses completed pipeline stages. Some incomplete stages must run again.
Assembly copies manifests and exported MTZs; it does not edit source results.

Default paths:

- Input: `~/.codex/pymda_data`
- Existing results: `~/.codex/pymda_output`
- Jobs and queue: `~/.codex/pymda_output/jobs`
- Programs: `/home/aa1/ccp47`

Override them with `PYMDA_INPUT_ROOT`, `PYMDA_OUTPUT_ROOT`,
`PYMDA_JOBS_ROOT`, and `PYMDA_CCP4` before starting. The process needs write
access to the jobs directory. Selectable directories remain within these roots.
Use `python run_web.py --port 8001` to change the port.

The API provides `GET /api/datasets`, `GET/POST /api/jobs`,
`GET /api/jobs/{id}`, and `POST /api/jobs/{id}/cancel` or `/resume`.
`GET /openapi.json` describes the request fields. Interactive API documentation
is disabled because it requires JavaScript. Files are available through
`GET /jobs/{id}/file?path=best/final.mtz`.

Merge selection compares finite crystal subsets, frame prefixes and resolution
cutoffs. It ranks eligible results by resolution, then overall merged I/σ(I),
then outer-shell CC½. Inspect completeness and the alternatives as well as the
selected result; this search does not implement the paper's SmRmerge algorithm
or prove an absolute optimum. Record whether scans represent independent
crystals or repeated exposures of one crystal.

Regression tests:

```bash
.pymda/bin/python -m unittest -v test_pymda test_web
```

The explicit real-data smoke test requires a running server, the default four
input scans, and existing crystal manifests. It submits a two-dataset assembly
and a 180-frame single-scan processing pilot; it can take several minutes and
retains output jobs:

```bash
.pymda/bin/python test_web_live.py
```
