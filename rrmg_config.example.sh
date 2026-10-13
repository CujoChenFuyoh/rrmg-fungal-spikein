#!/usr/bin/env bash
# RRMG site configuration.
#
# Copy this file to rrmg_config.sh, edit the absolute paths, and keep the
# resulting rrmg_config.sh private. It is ignored by Git.

# Directory containing the cloned RRMG code.
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Working directory containing input/, ref/, temp/, result/, and logs/.
WD="/absolute/path/to/SpikeIN"

# One integrated Conda environment containing all analysis programs.
CONDA_SH="${HOME}/miniconda3/etc/profile.d/conda.sh"
RRMG_CONDA_ENV="fungamr_meta"

# RandomReadsMG/BBMap and Java. Use absolute paths.
BBMAP_DIR="/absolute/path/to/bbmap"
RANDOMREADSMG_SH="${BBMAP_DIR}/randomreadsmg.sh"
JAVA_HOME="/absolute/path/to/java"
JAVA_BIN="${JAVA_HOME}/bin/java"

# Paths below may be absolute or relative to WD.
PLAN="input/rrmg/rrmg_permutation_plan.tsv"
BACKGROUND_PLAN="temp/rrmg_background/background_plan.tsv"
BACKGROUND_INDEX="temp/rrmg_background/background_index.tsv"
GENE_REF="ref/CYP51A_wt.fasta"
TARGET_GENOME="ref/Afum_Af293.fasta"

# Kraken2 is optional. Set RUN_KRAKEN=1 and provide a database to enable it.
RUN_KRAKEN=0
KRAKEN_DB="/absolute/path/to/kraken2_database"

# Runtime parameters. Keep these consistent with the PBS resource requests.
# Non-PBS runs default to one row at a time. Increase only after accounting for
# each row's thread and memory settings below.
LOCAL_JOBS=1

STEP1_THREADS=8
STEP1_BBMAP_XMX="32g"
BACKGROUND_SEED=1000
STEP1_FORCE=0

STEP2_KEEP_TARGET_READS=1
STEP2_FORCE=0

STEP3_THREADS=64
STEP3_SPADES_MEM_GB=450
STEP3_PHRED_OFFSET=33
STEP3_SKIP_QC=0
STEP3_SKIP_ASSEMBLY=0
STEP3_FORCE=0

STEP3B_THREADS=32
VAMB_MIN_CONTIG=2000
VAMB_MIN_BIN_SIZE=200000
VAMB_SEED=42
VAMB_COMMAND="vamb"
VAMB_USE_CUDA=0
VAMB_FORCE=0

STEP45_THREADS=8
STEP45_BIN_SIZE=10000
STEP45_MAKE_BINS=1
STEP45_FORCE=0
