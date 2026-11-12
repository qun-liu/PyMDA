# Processed-output discovery validation

Validated on November 11, 2026.

All 40 regression tests passed with `python3 -m unittest -v test_pymda test_web`.
Coverage includes parent output and job directories, crystal and export folders,
distinct stable IDs for repeated runs, one default version per scan, exact export
copying, unavailable results, moved exports, symlinked crystal directories, and
assembly defaults consistent with integration limits. Existing single-crystal,
wedge, bounded-subset and worker recovery tests also passed.

Live tests passed against the running localhost app using:

```bash
python3 test_web_live.py --source /home/aa1/.hermes/codex/test1 --input-dir /home/aa1/.codex_bak/pymda_data --skip-process --resolution 3.0 --output-dir /tmp/pymda-discovery-checks
python3 test_web_live.py --source /tmp/pymda-selected-output-k9sy3ni_ --input-dir /home/aa1/.codex_bak/pymda_data --skip-process --output-dir /tmp/pymda-discovery-checks
```

These verified parent-folder discovery, preferred result selection, single- and
two-dataset assembly, intensity comparison, compatibility-group selection,
selection provenance, final MTZ downloads, subset plans and wedge MTZ downloads.
Source manifests and exported MTZs retained their original SHA-256 hashes.
Existing integrations were reused; no raw-image reintegration was necessary.

The test1 scan assembled at 3.0 Å with 96.5% overall completeness and overall
I/sigma 17.4. The two-scan assembly reached 2.5 Å with 95.2% completeness and
overall I/sigma 19.7.

The user removed the former validation directory during verification. The app
was restarted with its persistent queue in `web_state/queue.sqlite3`, separate
from result directories. Five known completed jobs were recovered from their
remaining output artifacts, including the user's original processing job and
assembly 190ec855ca174f0393b4180b7fc3126c. That assembly's final MTZ download was
verified against the on-disk file byte-for-byte. Removed validation outputs
were not recreated.

The app is running at http://127.0.0.1:8000 with test1 as its initial processed
source and output root. Parent, crystals, individual crystal and prefix-folder
pages were checked through HTTP and returned successfully. No jobs were queued
or running at the end of verification. Launch settings are documented in WEB.md.

Live logs: `/tmp/pymda-live-test1-discovery.log` and
`/tmp/pymda-live-multi-discovery.log`.
