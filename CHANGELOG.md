# Changelog

## 0.3.0 — 2026-10-14

- BREAKING: rrmg_submit.sh executes rows directly; it no longer calls qsub or
  accepts qsub options. Submit user-owned PBS scripts explicitly with qsub.
- Both Bash entry points and rrmg run share a scheduler-free single-row dispatcher.
  Add step4 and step5 separately; keep step45 and optional step3b.
- Select original run_index/bg_index values, validating duplicates/missing IDs
  before execution. No synthetic PBS variables or PBS files in the core runner.
- Resolve all programs from the activated environment. Remove Step3 conda run;
  datasets and Java use the same environment. Explicit Java/BBMap overrides remain.
- Add manifest/download/gatekeeper instructions and local, ordinary PBS and
  PBS-array examples. Legacy PBS wrappers delegate to the same core.
- Conservative pilot resource defaults; bound Step2/test Java heap and threads.
  Allow Step2 to regenerate target outputs when rerunning after partial results.
- Add portable dispatch, environment selection and genome-preparation tests.

Read-count semantics, experimental seeds, result layouts and mutation-calling
limitations are unchanged. A successful collector is not a completeness audit.

## 0.2.0 — 2026-10-13

- Add PEP 517/621 packaging: `pip install .`, Git URL installs, source archives and wheels.
- Add `rrmg` entry point (`init`, `export`, `path`, `run`, `doctor`, `script`).
- Keep the original flat script layout as the single source of analysis code.
- `rrmg run` delegates to the existing local runner; no scheduler or new analysis engine is introduced.
- Support binding an existing checkout/config or exporting a versioned core snapshot without touching scientific inputs/results.
- Export explicit core files with SHA-256 manifest; exclude private config and sequencing data from distributions.
- Fix supplementary collector lookup when code and working directory differ.
- Add packaging/dispatch tests and Linux CI install checks.

Scientific simulation parameters, read-count semantics, analysis algorithms,
resume checks and existing force behavior are unchanged. This release does not
resolve previously observed read-count anomalies or add fungAMR variant calling.

## 0.1.0

- Original split workflow with configured PBS wrappers and a local runner.
