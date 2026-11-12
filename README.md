# PyMDA

Microcrystal diffraction data assembly using Python, DIALS and CCP4.
Original method: Takemaru et al., J. Appl. Cryst. (2020),
https://doi.org/10.1107/S160057671901673X.

This version uses Python 3 and current DIALS `.expt`/`.refl` files. It sources
CCP4's `bin/ccp4.setup-sh` automatically. The default installation is found with
`which aimless`; override it with `--ccp4` or `PYMDA_CCP4`. `--ccp4` accepts an
installation, its `bin` directory, or a parent containing `ccp4-*`.

```bash
python3 pymda.py --discover
python3 pymda.py --run_dials --run_mda --thread 2 --wedges 4
```

For this test collection, the finer cutoff search used:

```bash
python3 pymda.py --run_dials --run_mda --thread 2 --wedges 4 \
  --merge_workers 4 \
  --resolutions 2.0 2.2 2.22 2.24 2.25 2.26 2.27 2.28 2.3 2.4 2.5 3.0
```

`--merge_workers` controls independent subset jobs, separately from DIALS
image-processing workers. The search reuses completed stages and scaling trials.

Defaults: tools from the CCP4 installation found by `which aimless`, recursively discovered rotation master
files under `~/.codex/pymda_data`, output under `~/.codex/pymda_output`.
Broken links, duplicate master-file links and raster files are excluded. HDF5
master files must retain accessible companion data files. Raw inputs are never
modified. The controlling Python needs only the standard library; DIALS/CCP4
supply their own scientific runtime.

The pipeline imports, finds spots, indexes without imposing a space group,
refines and integrates each scan to 1.8 A. `--spg` and `--unit_cell` accept known
symmetry/cell constraints. Each crystal is exported at cumulative frame prefixes
(`--wedges 4`: 25%, 50%, 75%, 100%) without reintegrating the images. Prefixes
are assessed with Pointless/Aimless at a common cutoff (`--reso_cchalf 2.5`);
the prefix with strongest overall merged I/sigma(I), provided its overall
CC-half is at least 0.3 when possible, is selected for subset comparisons.
This is a modern search inspired by PyMDA, rather than a literal implementation
of the paper's SmRmerge frame-rejection procedure. The legacy interactive
rejection flags have been replaced by an automatic, auditable search.

`--split_wedges 3` additionally divides a 1,800-frame scan into consecutive
wedges 1-600, 601-1200 and 1201-1800. Each wedge receives its own intensity MTZ
and Pointless/Aimless assessment. The default, 1, keeps the scan unsplit.
If the frame count is not divisible by the wedge count, the extra frames go
to the earliest wedges; no frame appears in two wedges. Boundaries use the
original scan or pilot range, with DIALS-trimmed endpoints recorded separately.
Completed integrations can be reused when changing this export setting.

For collections of up to ten datasets, all subsets of two or more crystals are scaled across `--resolutions`
using I/sigma-selected prefixes, CC-half-selected prefixes (as in the paper),
and full-scan prefixes; identical trials are computed once. The cutoffs default to 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.0 A. Pointless checks symmetry
and indexing compatibility, rather than copying an assumed indexing. Aimless
uses its standard error model: I/sigma(I) is not inflated by disabling error
corrections. Quality thresholds default to 90% overall completeness (following the paper), 80% outer-shell completeness, outer-shell
CC-half >= 0.3 and outer-shell merged I/sigma(I) >= 1.0. Eligible candidates are
ranked by highest resolution (smallest A), then overall I/sigma(I), then outer
CC-half. Resolution, signal and completeness are separate tradeoffs: this finite
search cannot establish an absolute optimum. Inspect shell statistics and
anisotropy before structural interpretation. Free-R flags are created once for
the chosen output and should be retained during later refinement.

`best/final.mtz` contains truncated amplitudes and 5% Free-R flags;
`best/scaled.mtz` retains merged intensities and anomalous columns.
`alternatives/` preserves the strongest-signal merge and, when available, the
highest-resolution merge reaching 95% overall completeness. An
`extended_resolution` alternative preserves the measured high-resolution data
with supported signal even when completeness fails the selection thresholds;
its statistics must be considered when using it. `pareto.json` and
`alternatives/pareto/` retain nondominated resolution/signal choices: improving
resolution or overall I/sigma further requires sacrificing the other.
`--crystal_identity` records known repeat/independent crystal provenance; its
default is unknown. Unit-cell agreement alone does not establish crystal identity.
`resolution_tradeoff.png` plots the resolution/signal/completeness tradeoff.
`best.json`, `trials.json`, `trials.csv`, individual `statistics.json`, tool command JSON and
console logs record the selection and provenance. Exact repeats retain the
existing Free-R flags. A changed reflection dataset creates a new final MTZ;
retain previously assigned flags if continuing an existing refinement project. `--nres` optionally supplies
the known number of residues in the asymmetric unit to ctruncate.

Completed stages, crystal manifests and successful scaling trials are reusable.
`--exports_only` recovers exports from completed integrations, using their actual
scan bounds (DIALS may trim endpoints), without rereading raw images. To
repeat assembly alone:

```bash
python3 pymda.py --run_mda --resolutions 2.0 2.1 2.2 2.3 2.5
```

Use a new output directory when changing integration/prefix settings. `--frames`
can limit a pilot run; it changes the scientific sampling and must be documented.
Individual crystal/trial failures are logged; a run exits unsuccessfully if no
result meets the stated quality thresholds. `--subset_search auto` accepts one
or many datasets: up to ten it compares all subsets; above ten it uses a bounded
search with at most `--max_subsets 64` candidates per prefix policy. Candidates
include the full collection, available compatibility groups, leave-one-out
trials for up to sixteen individually weakest scans, and growing subsets ranked
by individual signal with CC-half support. `search.json` records every candidate
and the actual method. This bounded search does not explore every combination.
Use `--subset_search bounded --max_subsets N` to change the budget, or
`--subset_search exhaustive` for a collection of at most ten datasets.
Standalone `uc_rejection.ucr` and `picc_rejection.pr` clustering
utilities remain available; the CLI subset search does not depend on SciPy.

To finalize a successfully scaled result even when no trial passes the quality
thresholds, add `--skip_quality_checks`. Measured statistics and the original
`eligible` pass/fail value remain in `trials.json` and `best.json`; the report
and `best.json` record that quality checks were bypassed. Actual scaling or
finalization failures still stop processing.

For one merge of all selected datasets, without searching alternative subsets,
prefixes or resolution cutoffs, use one cutoff and a bounded search of one:

```bash
python3 pymda.py --run_mda --dir_output /path/to/processed/output \
  --dataset_ids ID1 ID2 --resolutions 3.0 --prefix_policy full \
  --subset_search bounded --max_subsets 1 --skip_quality_checks
```

Scaling still runs and writes `trials.json`. If that file is empty, inspect the
trial `failure.json` files and Pointless/Aimless logs; threshold bypass cannot
recover a failed scaling operation.

```bash
python3 -m unittest -v test_pymda
```

## Local web application

Run `python3 run_web.py` and open http://127.0.0.1:8000.
The FastAPI interface uses Python-generated pages with no Jinja2 or JavaScript,
and provides dataset selection, queued processing/assembly, cancellation,
recovery, logs, statistics, plots and downloads. No authentication is required;
access is restricted to localhost. See [WEB.md](WEB.md) for configuration,
API routes and regression/real-data tests.

Single datasets can now be finalized with `--run_dials --run_mda` and a single
selected master file. To compare exported datasets, use `--compare` against
an output directory containing crystal manifests; optional screening controls
are `--cell_length_tolerance`, `--cell_angle_tolerance`, `--min_pair_cc`,
and `--min_common`, with `--reso_cchalf` setting the common comparison cutoff.
The web interface provides comparison, group selection and directory buttons;
see [WEB.md](WEB.md).

Use `--dataset_ids ID1 ID2` to restrict comparison or assembly to processed
manifest directory IDs. `selection.json` records included and excluded IDs.

The same workflow handles other rotation HDF5 collections using `--hdf5_path`
and a collection-specific `--dir_output`. Rerun the same command to resume:
completed imports, integrations, exports and successful scaling trials are
reused, and failed stages are retried. `--exports_only` also supports completed
pilot integrations, after verifying their original scan bounds and settings.
