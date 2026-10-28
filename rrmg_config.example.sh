#!/usr/bin/env bash
# RRMG site configuration.
#
# Copy this file to rrmg_config.sh, edit the absolute paths, and keep the
# resulting rrmg_config.sh private. It is ignored by Git.

# Directory containing the cloned RRMG code.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Working directory containing input/, ref/, temp/, result/, and logs/.
WD="/absolute/path/to/SpikeIN"

# Activate fungamr_meta BEFORE running this workflow. Tools default to its PATH;
# explicit external tool overrides are supported. No second environment needed.
# CONDA_SH / RRMG_CONDA_ENV are optional conveniences for PBS examples only.
CONDA_SH=""
RRMG_CONDA_ENV="fungamr_meta"

# Optional tool overrides. Normally leave empty (auto-resolve from active PATH).
# A separately installed BBMap can be selected without changing environments.
# Check both launchers, not just the environment name; see docs/BBMAP.zh-CN.md.
# Explicit launcher overrides win over BBMAP_DIR/PATH; clear old overrides
# before switching back to a Conda-managed BBMap installation.
BBMAP_DIR="${BBMAP_DIR:-}"
RANDOMREADSMG_SH="${RANDOMREADSMG_SH:-}"
RANDOMREADS_SH="${RANDOMREADS_SH:-}"
JAVA_HOME="${JAVA_HOME:-}"
JAVA_BIN="${JAVA_BIN:-}"

# Paths below may be absolute or relative to WD.
PLAN="input/rrmg/rrmg_permutation_plan.tsv"
BACKGROUND_PLAN="temp/rrmg_background/background_plan.tsv"
BACKGROUND_INDEX="temp/rrmg_background/background_index.tsv"
GENE_REF="ref/CYP51A_wt.fasta"
TARGET_GENOME="ref/Afum_Af293.fasta"

# Kraken2 is optional. Set RUN_KRAKEN=1 and provide a database to enable it.
RUN_KRAKEN=0
KRAKEN_DB="/absolute/path/to/kraken2_database"

# Conservative fallback defaults, used only when the command omits -t/--threads
# or (for Steps 1-3) -m/--memory. They are NOT guaranteed sufficient for a full
# metagenome. Command-line memory has no artificial maximum, but it does not
# request scheduler resources; match the PBS/Slurm allocation and leave overhead.
# Non-PBS runs default to one row at a time. Increase only after accounting for
# each row's thread and memory settings below.
LOCAL_JOBS=1

STEP1_THREADS=4
STEP1_BBMAP_XMX="4g"
BACKGROUND_SEED=1000
STEP1_FORCE=0

# One thread keeps target simulation deterministic with a fixed BBMap seed.
STEP2_THREADS=1
STEP2_BBMAP_XMX="4g"
STEP2_KEEP_TARGET_READS=1
STEP2_FORCE=0

STEP3_THREADS=4
STEP3_SPADES_MEM_GB=16
STEP3_PHRED_OFFSET=33
# 1: pass --only-assembler (skip BayesHammer, but still run fastp QC).
# 0: enable native metaSPAdes read error correction before assembly.
# Keep this choice consistent across compared runs. A mode change requires
# rebuilding Step3 with --force; existing assemblies are never silently relabelled.
STEP3_ONLY_ASSEMBLER=1
STEP3_SKIP_QC=0
STEP3_SKIP_ASSEMBLY=0
STEP3_FORCE=0

STEP3B_THREADS=4
VAMB_MIN_CONTIG=2000
VAMB_MIN_BIN_SIZE=200000
VAMB_SEED=42
VAMB_COMMAND="vamb"
VAMB_USE_CUDA=0
VAMB_FORCE=0

STEP45_THREADS=4
# Step4 maps both mixed metagenome reads and retained target-only reads to GENE_REF.
STEP4_MAP_TARGET_ONLY=1
STEP45_BIN_SIZE=10000
STEP45_MAKE_BINS=1
STEP45_FORCE=0
