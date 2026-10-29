# PyMDA

Microcrystal diffraction data assembly using Python, DIALS and CCP4.
Original method: Takemaru et al., J. Appl. Cryst. (2020),
https://doi.org/10.1107/S160057671901673X.

This version uses Python 3 and current DIALS `.expt`/`.refl` files. It sources
CCP4's `bin/ccp4.setup-sh` automatically; `--ccp4` accepts either an installation
or a parent containing `ccp4-*`. No global PATH edits are needed.

```bash
python3 pymda.py --discover
python3 pymda.py --run_dials --run_mda --thread 2 --wedges 4
```
python3 pymda.py --run_dials --run_mda --thread 2 --wedges 4 \
  --merge_workers 4 \
  --resolutions 2.0 2.2 2.3 2.4 2.5 3.0
```

`--merge_workers` controls independent subset jobs, separately from DIALS
image-processing workers. The search reuses completed stages and scaling trials.

Defaults: tools under `/home/aa1/ccp47`, recursively discovered rotation master
files under `~/pymda_data`, output under `~/pymda_output`.
Broken links, duplicate master-file links and raster files are excluded. HDF5
master files must retain accessible companion data files. 

The pipeline imports, finds spots, indexes refines and integrates. `--spg` and `--unit_cell` accept known symmetry/cell constraints. Each crystal is exported at cumulative frame prefixes
(`--wedges 4`: 25%, 50%, 75%, 100%) without reintegrating the images. Prefixes
are assessed with Pointless/Aimless at a common cutoff (`--reso_cchalf 2.5`);
the prefix with strongest overall merged I/sigma(I), provided its overall
CC-half is at least 0.3 when possible, is selected for subset comparisons.

All subsets of two or more crystals are scaled across `--resolutions`
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
can limit a pilot run; it changes the scientific sampling.
Individual crystal/trial failures are logged; a run exits unsuccessfully if no
supported multi-crystal result exists. Exhaustive subset comparisons are limited
to ten crystals. Standalone `uc_rejection.ucr` and `picc_rejection.pr` clustering
utilities remain available; the CLI subset search does not depend on SciPy.

```bash
python3 -m unittest -v test_pymda
```

## Local web application
python3 -m venv ~/.pymda
source ~/.pymda/bin/activate

#install FastAPI using pip, preferably using an isolated environment.

pip install "fastapi[standard]"

Run `~/.pymda/bin/python run_web.py` and open http://127.0.0.1:8000.
The FastAPI interface uses Python-generated pages 
and provides dataset selection, queued processing/assembly, cancellation,
recovery, logs, statistics, plots and downloads. No authentication is required;
access is restricted to localhost. See [WEB.md](WEB.md) for configuration,
API routes and regression/real-data tests.
