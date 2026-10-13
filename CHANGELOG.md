# Changelog

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
