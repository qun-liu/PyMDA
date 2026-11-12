# Local PyMDA interface

Start from the repository directory:

```bash
python3 run_web.py
```

Open http://127.0.0.1:8000. Stop the server with Ctrl-C. The interface uses
FastAPI and Python-generated HTML/CSS. Directory buttons, forms and five-second page refreshes
work without JavaScript or Jinja2. There is no authentication; the server binds
to loopback and accepts local requests only.

The workspace recursively discovers diffraction master HDF5 files and checks
frame counts and companion files with CCP4 Python. Select scans to process,
process and assemble, or select existing processed datasets for assembly alone.
Advanced parameters include pilot frame limits, geometry, resolution cutoffs,
quality thresholds and crystal provenance. Use a pilot before a full scan run.

Under **Processing parameters**, set **Consecutive wedges per scan** to the
desired count (default 1). For a 1,800-frame single-crystal dataset, 3 creates
wedges 1-600, 601-1200 and 1201-1800. Each wedge has a downloadable intensity
MTZ and quality statistics in the job's **Consecutive scan wedges** table.
The count applies separately to each selected scan and to the limited range
when a pilot frame limit is supplied. Unequal divisions differ by at most one
frame, with extra frames assigned to the earliest wedges. Requested bounds
and the frames retained by DIALS are both shown when endpoints are trimmed.
**Cumulative frame prefixes** remains the separate prefix-assessment setting.

One or many datasets can be processed, compared or assembled. **Dataset subset
search** defaults to automatic: exhaustive through ten datasets and bounded for
larger collections. **Maximum bounded-search subsets** controls the candidate
budget (default 64 per prefix policy). Bounded candidates include the full
collection, available compatibility groups, removal of individually weak scans,
and progressively larger subsets ranked by individual signal. Results link to
`search.json` so the evaluated choices are visible. For other collections,
choose their input and output directories; use a separate output location for
each collection. Failed or interrupted jobs can be resumed with their original
parameters, reusing completed stages and retrying failures.

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
access to the jobs directory. These paths are starting locations, not access
limits. Input and processed-result folders can be anywhere the workstation
account has read and traversal access. Output folders can be anywhere it has
write and traversal access; missing output folders are created when a valid job
is submitted. Directory browsing and workspace visits do not create folders.
Deleted or inaccessible saved/default folders fall back to an accessible
starting location. Explicitly selected input and processed-result folders must
exist.
Use `python run_web.py --port 8001` to change the port.

The API provides `GET /api/datasets`, `GET /api/processed-datasets?source_dir=...`, `GET/POST /api/jobs`,
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
python3 -m unittest -v test_pymda test_web
```

The explicit real-data smoke test requires a running server, readable input
scans, and existing processed results. It verifies discovery, comparison,
single-dataset assembly, and two-dataset assembly when two completed scans are
available. It also runs a 180-frame processing pilot unless `--skip-process`
is supplied; it can take several minutes and retains output jobs:

```bash
python3 test_web_live.py --source /path/to/output --input-dir /path/to/raw/data
```

Use `--output-dir /path/to/test/output` to keep the new test jobs separate and
`--resolution 3.0` for data integrated to 3 Å.

Single-dataset finalization is available through **Process and assemble** or
**Scale / assemble selected datasets** with one dataset selected. Processing
alone retains exported prefixes for later comparison or finalization.

Use **Compare selected datasets** on processed data to calculate Niggli-reduced
cell differences and Pearson intensity correlations after Pointless aligns
indexing against a reference. The comparison uses merged, non-anomalous common
reflections at the specified cutoff. Defaults are 2% cell-length difference,
2 degrees angle difference, correlation >= 0.9, and at least 100 common
reflections. These are adjustable screening criteria, not guarantees of
isomorphism. Reduced-cell comparisons can be sensitive near reduction boundaries.
The heatmap reports missing correlations and overlap counts. Suggested groups
require every pair to pass; select a group, adjust individual checkboxes, then
submit assembly. Comparison and selection JSON preserve the evidence and choices.

Use **Browse data directory**, **Browse processed results**, or **Choose output
directory** to select accessible folders anywhere on the server workstation.
The browser offers parent-folder,
home-folder and filesystem-root navigation, a folder-path field, and **Create
folder**. Click **Use this directory** to apply the selection to processing and
assembly. The app remembers it in this browser for subsequent workspace visits.
Results saved outside the default output root can also be selected for assembly.
Outputs use `<chosen directory>/<job id>/output`;
the durable queue stays under `PYMDA_JOBS_ROOT`. The output browser can create
folders, or you can enter a new output path directly in the job form.

Click **Find processed datasets in this output folder** beside the output
destination to list the results already saved there. The search finds
`crystals/<dataset>/manifest.json` inside nested job folders, so selecting the
parent output directory is sufficient. Job folders, their `output` directories,
`crystals` directories, individual crystal directories and their export folders
also work. Each row shows the result location, integration limit and readiness;
**Open this result directory** opens the scans from that run.

Separate runs of the same scan remain individually selectable. One complete
version per raw master file is checked by default, preferring full scans over
pilots, more frames, a finer integration limit, then the newest result. Incomplete
manifests or missing exports appear as unavailable while complete results remain
selectable. Assembly and comparison defaults respect the selected results'
integration limits; check them again if you change the dataset checkboxes.
The app remembers the processed-data directory in this browser. Moved results
can be reused when their exported MTZs remain in their crystal directories.
API clients should submit catalog `id` values as `dataset_ids`; these distinguish
separate runs and remain stable when additional jobs are added to the same folder.

For this checkout and supplied data, launch with:

```bash
PYMDA_INPUT_ROOT=/home/aa1 \
PYMDA_START_INPUT_DIR=/home/aa1/.codex_bak/pymda_data \
PYMDA_START_SOURCE_DIR=/home/aa1/.hermes/codex/test1 \
PYMDA_OUTPUT_ROOT=/home/aa1/.hermes/codex/test1 \
PYMDA_JOBS_ROOT="$PWD/web_state" python3 run_web.py
```

Then browse to `/home/aa1/.codex_bak/pymda_data` for input. Hidden folders can
also be entered in the directory text field. Directory choices persist when
rediscovering scans or processed datasets.
Keep `web_state` to retain the job queue and history; result directories can be
selected independently of that state folder.

Scientific correlation regression checks run in the CCP4 runtime:

```bash
PYMDA_TEST_MTZ=/absolute/path/to/exported/integrated.mtz \
/home/aa1/ccp47/ccp4-9/bin/ccp4-python -m unittest -v test_compare_scientific
```

This verifies common-intensity correlation, three-dataset complete-link grouping,
and unavailable correlations when the required overlap exceeds measured data.
