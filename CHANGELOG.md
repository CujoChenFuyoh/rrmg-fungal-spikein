# Changelog

## 0.3.4 — 2026-10-30

- Add independently resumable contig cs/cg PAF and raw/Q20 base counts for both
  metagenomic and synthetic target-only reads. Step3/4 generate these by default;
  existing projects can backfill them from contigs/BAMs without rebuilding reads,
  QC, assembly or existing mapping results.
- Validate reference, generation, source identity, options, tool versions and
  output checksums for the new derived evidence. Missing legacy generation stays
  UNVERIFIED; zero coverage is explicit and missing BAMs are not biological zeros.
- Collect portable evidence, per-contig nucleotide differences, full/fragmented
  recovery and review-only discordance flags. Add --require-base-evidence for a
  strict evidence-complete snapshot. Do not infer false positives or FungAMR calls
  solely from alignment breadth. Existing R coverage plotting remains compatible.

## 0.3.3 — 2026-10-28

- Extend Step4 to map both mixed metagenome reads and retained target-only
  reads to the target gene, producing independent BAM, PAF, per-base depth and
  coverage summaries. Resume each read set independently so existing audited
  mixed-read outputs are not recomputed merely to add target-only results.
- Replace post-analysis remapping with a read-only collector that validates and
  packages existing Step3-5 outputs into one portable `figure_data/` snapshot.
- Simplify the runbook ending to two stages: collect/audit figure inputs on
  Linux, then download the bundle and render both R figure families on Windows.
- Add a Windows PowerShell plotting launcher and retain fixed, separate
  `figure_data/` input and `figure_output/` rendering directories.
- Collect mixed-read and target-only-read Step4 summaries separately, then
  require exact plan membership before creating a complete figure-data bundle.

## 0.3.2 — 2026-10-19

- Consolidate user documentation into quick start, installation, workflow,
  figure and troubleshooting references; remove personal path defaults and
  unversioned manual validation claims from user instructions.
- Keep PAF generation enabled for complete read/contig alignment figures.
  Mapping reuse and unified reference/parameter fingerprints remain pending.
- Use conservative single-row pilot resources without implicit force rebuilds.
- Remove the five fixed-stage root PBS array wrappers and the obsolete
  target-only figure patch. Keep the two generic templates in `pbs_example/`
  and the unified alignment-input preparation, with default PAF generation.
- Rebuild generated wheel resources from the current core asset list, so an
  older local build cache cannot reintroduce removed workflow entries.
- Count reference sequence bases without FASTA whitespace in Step3; LF and
  CRLF inputs now give the same contig breadth. Existing summaries must be
  recalculated; reads, assembly, mapping thresholds and PAF spans are unchanged.
- Reject duplicate seeds/designs/run IDs and the inconsistent `copy_number`
  mode before execution. Keep supported `read_fraction` calculations unchanged.
- Resolve the optional VAMB log environment without requiring an exported
  `RRMG_CONDA_ENV`; detect/configure Windows R without a personal library default.
- Preserve existing site deployment files, plans and results. No Git history
  rewrite, license selection, release tag or automatic push is performed.

## 0.3.1 — 2026-10-15

- Document Conda-managed versus external BBMap, both required launchers,
  fixed-version upstream archive downloads, override precedence and PBS checks.
  Include the installation guide in pip/export assets; no environment is
  installed or changed automatically.
- Default Step3 to `STEP3_ONLY_ASSEMBLER=1` (metaSPAdes `--only-assembler`)
  while keeping fastp QC. Set `0` to enable native BayesHammer. Record assembly
  mode and refuse reuse of mismatched/unrecorded modes without an explicit
  rebuild. Document method comparability and re-collection in the runbook.
- Normalize the locally CRLF-encoded runbook back to LF for Bash portability.
- Let each invocation override per-row threads with `-t/--threads`. Add
  `-m/--memory GB` for the BBMap Java heap in Steps 1/2 and metaSPAdes in
  Step 3, with no workflow-imposed maximum; scheduler memory remains external.
- Bundle the step-by-step local runbook, R figures, WSL-to-Windows R launcher
  and alignment-input preparation in both the repository and pip/export assets.
- Use CODE for all scripts, with fixed figure_data/figure_output paths under WD.
  Render/map in staging; replace only verified tool-owned output directories
  after success, retaining previous outputs on failure and rejecting user edits.
- Validate empty/missing annotation paths before mapping; document the three
  gene PAF inputs and the separate target-only gene versus genome mapping.
- Make depth overlays explicit: pale grey target-only fill, stronger orange
  mixed-QC foreground, dashed/solid outlines, covered-base counts and separate
  read-track opacity. Test identical/partial overlap and true-zero profiles.
- Pin repository text files, including R scripts and package metadata, to LF
  on checkout even when Windows Git has core.autocrlf=true. Exclude figure
  outputs from Git and add a source line-ending regression check.
- Replace byte-concatenated background/target gzip streams with deterministic
  paired interleaving into one gzip member. Record exact background, target and
  total counts plus output sizes in each run's read_mix_audit.tsv.
- Use BBMap RandomReads for an exact whole-target pair total; uniformly
  subsample oversized RandomReadsMG backgrounds to the designed count.
- Fail closed when fastp did not consume exactly twice the planned pair count,
  when target-only pairs differ from design, or when rebuilt reads would be
  combined with stale QC, assembly or mapping outputs. Legacy unaudited mixes
  must be rebuilt rather than silently reused.
- Assign a new mix generation ID on every Step2 rebuild and propagate it through
  QC, assembly and Step4/5 provenance markers, so file touching or preserved
  mtimes cannot make stale downstream outputs look current.

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
