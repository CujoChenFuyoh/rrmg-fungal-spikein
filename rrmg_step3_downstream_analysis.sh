#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_MIXER="${FASTQ_MIXER:-$SCRIPT_DIR/rrmg_mix_paired_fastq.py}"

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step3_downstream_analysis.sh plan.tsv ROW_INDEX -w /path/to/project \
    -gene ref/CYP51A_wt.fasta --db /path/to/kraken_db -t 64 -mem 450

Purpose:
  Run downstream analysis from already-merged reads:
    optional fastp -> optional kraken2 -> optional metaSPAdes -> minimap2 to gene/species refs.

Interpretation:
  This split Step 3 measures target-gene detection, alignment identity, and
  breadth. It currently creates an empty VCF placeholder and reports VarN=0;
  it does NOT call or annotate CYP51A/fungAMR mutations.

Options:
  -gene FILE             Target gene FASTA for contig mapping
  --species-ref FILE     Optional target species/genome FASTA for contig mapping
  --db DIR               Kraken2 DB
  -species NAME          Species name to extract from Kraken2 report; default from plan
  -t INT                 Threads, default 24
  -mem INT               metaSPAdes memory GB, default 250
  --only-assembler       Skip BayesHammer read correction (default); keep fastp QC
  --with-error-correction  Enable native BayesHammer before assembly
  --skip-qc
  --skip-kraken
  --skip-assemble       Reuse an existing non-empty asm/contigs.fasta
  --force               Replace existing QC/Kraken/assembly outputs
  --phred-offset INT     metaSPAdes PHRED offset, default 33
USAGE
}

PLAN="${1:-}"
ROW_INDEX="${2:-${PBS_ARRAY_INDEX:-}}"
if [[ -z "$PLAN" || -z "$ROW_INDEX" ]]; then usage; exit 1; fi
shift 2 || true

WD=""
GENE_REF=""
SPECIES_REF=""
KRAKEN_DB=""
SPECIES_NAME=""
THREADS=24
SPADES_MEM_GB=250
ONLY_ASSEMBLER=1
SKIP_QC=0
SKIP_KRAKEN=0
SKIP_ASSEMBLE=0
FORCE=0
PHRED_OFFSET=33

# Resolve analysis programs from the caller's active environment/PATH.

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WD="$2"; shift 2;;
    -gene) GENE_REF="$2"; shift 2;;
    --species-ref) SPECIES_REF="$2"; shift 2;;
    --db) KRAKEN_DB="$2"; shift 2;;
    -species) SPECIES_NAME="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -mem) SPADES_MEM_GB="$2"; shift 2;;
    --only-assembler) ONLY_ASSEMBLER=1; shift 1;;
    --with-error-correction) ONLY_ASSEMBLER=0; shift 1;;
    --skip-qc) SKIP_QC=1; shift 1;;
    --skip-kraken) SKIP_KRAKEN=1; shift 1;;
    --skip-assemble) SKIP_ASSEMBLE=1; shift 1;;
    --force) FORCE=1; shift 1;;
    --phred-offset) PHRED_OFFSET="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Missing -w}"
WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"
[[ -s "$PLAN" ]] || { echo "[ERROR] Plan not found: $PLAN" >&2; exit 1; }
[[ -s "$FASTQ_MIXER" ]] || { echo "[ERROR] FASTQ audit utility not found: $FASTQ_MIXER" >&2; exit 1; }
if [[ -z "$KRAKEN_DB" ]]; then SKIP_KRAKEN=1; fi

merge_intervals_cov_bp() {
  sort -k1,1n -k2,2n | awk 'BEGIN{cov=0;s=-1;e=-1} NR==1{s=$1;e=$2;next} {if($1<=e){if($2>e)e=$2}else{cov+=(e-s);s=$1;e=$2}} END{if(s!=-1)cov+=(e-s);print cov}'
}

marker_matches() {
  local marker="$1" expected="$2"
  [[ -s "$marker" && "$(tr -d '\r\n' < "$marker")" == "$expected" ]]
}

write_marker() {
  local marker="$1" value="$2" temporary="${1}.tmp.$$"
  printf '%s\n' "$value" > "$temporary"
  mv -f -- "$temporary" "$marker"
}

kraken2_species_abundance() {
  local rpt="$1" sp="$2"
  if [[ -z "$sp" || ! -s "$rpt" ]]; then echo -e "NA\tNA\tNA\tNA"; return 0; fi
  awk -v sp="$sp" '
    function norm(s,t){t=tolower(s);gsub(/_+/," ",t);gsub(/[[:space:]]+/," ",t);sub(/^ +/,"",t);sub(/ +$/, "", t);return t}
    BEGIN{spn=norm(sp);bestS="";bestS_clade=-1;bestAny="";bestAny_clade=-1}
    {name="";for(i=6;i<=NF;i++){name=(name==""?$i:name" "$i)};sub(/^[[:space:]]+/,"",name);nn=norm(name);if(index(nn,spn)>0){clade=$2+0;if($4=="S"){if(clade>bestS_clade){bestS_clade=clade;bestS=$0}}else{if(clade>bestAny_clade){bestAny_clade=clade;bestAny=$0}}}}
    END{line=(bestS!=""?bestS:bestAny);if(line==""){print "0.000000\t0\t0\tNA";exit};n=split(line,a,"\t");if(n<5)n=split(line,a,/[[:space:]]+/);printf "%.6f\t%s\t%s\t%s\n",a[1],a[2],a[3],a[5]}' "$rpt"
}

run_fastp() {
  local in1="$1" in2="$2" out1="$3" out2="$4" js="$5" ht="$6"
  mkdir -p "$(dirname "$out1")" "$(dirname "$out2")"
  fastp -w "$THREADS" -i "$in1" -I "$in2" -o "$out1" -O "$out2" -j "$js" -h "$ht" 2> "${out1}.fastp.log"
  [[ -s "$out1" && -s "$out2" ]] || { echo "[ERROR] fastp failed: ${out1}.fastp.log" >&2; exit 1; }
}

run_kraken2() {
  local in1="$1" in2="$2" rpt="$3" out="$4"
  mkdir -p "$(dirname "$rpt")" "$(dirname "$out")"
  kraken2 --db "$KRAKEN_DB" --paired "$in1" "$in2" --threads "$THREADS" --use-names --report-zero-counts --memory-mapping --report "$rpt" --output "$out" > "${out}.kraken2.log" 2>&1
  [[ -s "$rpt" ]] || { echo "[ERROR] kraken2 failed: ${out}.kraken2.log" >&2; exit 1; }
}

run_metaspades() {
  local in1="$1" in2="$2" outdir="$3"
  local -a command=(metaspades.py -t "$THREADS" -m "$SPADES_MEM_GB" --phred-offset "$PHRED_OFFSET")
  [[ "$ONLY_ASSEMBLER" -eq 0 ]] || command+=(--only-assembler)
  command+=(-1 "$in1" -2 "$in2" -o "$outdir")
  mkdir -p "$outdir"
  printf '[SPADES]'; printf ' %q' "${command[@]}"; printf '\n'
  "${command[@]}" 2> "${outdir}/metaspades.log"
  [[ -s "${outdir}/contigs.fasta" ]] || { echo "[ERROR] metaSPAdes failed: ${outdir}/metaspades.log" >&2; exit 1; }
}

contigs_to_target_summary() {
  local sid="$1" contigs="$2" target="$3" ttype="$4" out="$5"
  if [[ ! -s "$contigs" || ! -s "$target" ]]; then echo -e "${sid}\t${ttype}\t${target}\t${contigs}\tN\t0\t0\t0\t0\tNA\t0" >> "$out"; return 0; fi
  local paf="${contigs}.to_${ttype}.paf" vcf="${contigs}.to_${ttype}.vcf"
  minimap2 -t "$THREADS" -x asm5 "$target" "$contigs" > "$paf" 2> "${paf}.minimap2.log"
  if [[ ! -s "$paf" ]]; then echo -e "${sid}\t${ttype}\t${target}\t${contigs}\tN\t0\t0\t0\t0\t${vcf}\t0" >> "$out"; return 0; fi
  local best_line best_alnlen best_ident target_len covbp breadth varn
  best_line=$(awk 'BEGIN{best=0} $12>=20{if($11>best){best=$11;line=$0}} END{print line}' "$paf")
  if [[ -z "$best_line" ]]; then echo -e "${sid}\t${ttype}\t${target}\t${contigs}\tN\t0\t0\t0\t0\t${vcf}\t0" >> "$out"; return 0; fi
  best_alnlen=$(echo "$best_line" | awk '{print $11}')
  best_ident=$(echo "$best_line" | awk '{if($11>0)printf "%.6f",100*$10/$11;else print "0"}')
  # FASTA formatting is not sequence: ignore whitespace, including CRLF's CR.
  target_len=$(awk 'BEGIN{len=0}/^>/{next}{gsub(/[[:space:]]/, "");len+=length($0)}END{print len}' "$target")
  covbp=$(awk '$12>=20{s=$8;e=$9;if(s>e){t=s;s=e;e=t};print s"\t"e}' "$paf" | merge_intervals_cov_bp)
  breadth=$(awk -v c="$covbp" -v L="$target_len" 'BEGIN{if(L>0)printf "%.6f",100*c/L;else print "0"}')
  : > "$vcf" || true
  varn=0
  echo -e "${sid}\t${ttype}\t${target}\t${contigs}\tY\t${best_alnlen}\t${best_ident}\t${breadth}\t${covbp}\t${vcf}\t${varn}" >> "$out"
}

ROW_ENV="$(mktemp "${WD}/temp/rrmg_downstream_row.XXXXXX")"
python3 - "$PLAN" "$ROW_INDEX" "$WD" > "$ROW_ENV" <<'PY'
import csv, re, shlex, sys
from pathlib import Path
plan=Path(sys.argv[1]); row_index=int(sys.argv[2]); wd=Path(sys.argv[3])
def safe(s):
    s=re.sub(r"[^A-Za-z0-9_.-]+","_",str(s).strip()); return s.strip('_') or 'NA'
def q(v): return shlex.quote(str(v))
with plan.open(newline='') as fh:
    chosen=None
    for row in csv.DictReader(fh, delimiter='\t'):
        if int(row.get('run_index','0') or 0)==row_index: chosen=row; break
if chosen is None: raise SystemExit(f'[ERROR] Row not found: {row_index}')
sample=chosen.get('sample') or chosen.get('run_id') or f'row{row_index}'
seed=chosen.get('seed') or row_index
run_id=chosen.get('run_id') or f'{sample}__Seed{seed}'
run_dir=wd/'temp'/'rrmg_permutation'/safe(run_id)
readlen=chosen.get('readlen') or '150'
species=chosen.get('species_name_for_kraken') or chosen.get('target_name') or ''
total_pairs=chosen.get('total_pairs') or chosen.get('computed_total_pairs') or ''
target_abundance=chosen.get('target_abundance','') or '0'
target_pairs=int(round(float(total_pairs)*float(target_abundance))) if total_pairs else ''
background_pairs=int(round(float(total_pairs)))-target_pairs if total_pairs else ''
vals={'RUN_ID':run_id,'SAMPLE':sample,'SEED':seed,'RUN_DIR':run_dir,'MIX_R1':run_dir/'reads'/f'{safe(sample)}.Seed{seed}_1.fastq.gz','MIX_R2':run_dir/'reads'/f'{safe(sample)}.Seed{seed}_2.fastq.gz','SPECIES_FROM_PLAN':species,'TARGET_ABUNDANCE':target_abundance,'TOTAL_PAIRS':total_pairs,'TARGET_PAIRS_EXPECTED':target_pairs,'BACKGROUND_PAIRS_EXPECTED':background_pairs,'READLEN':readlen}
for k,v in vals.items(): print(f'{k}={q(v)}')
PY
# shellcheck disable=SC1090
source "$ROW_ENV"; rm -f "$ROW_ENV"

if [[ -z "$SPECIES_NAME" ]]; then SPECIES_NAME="$SPECIES_FROM_PLAN"; fi
CONTIGS="${RUN_DIR}/asm/contigs.fasta"
MAP_SUM="${RUN_DIR}/map/target_call_summary.tsv"
SUCCESS_MARKER="${RUN_DIR}/STEP3_SUCCESS"
MASTER="result/rrmg_split_pipeline_summary.tsv"
MIX_AUDIT="${RUN_DIR}/design/read_mix_audit.tsv"
FASTP_JSON="${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}.fastp.json"
QC_PROVENANCE="${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}.mix_generation.txt"
KRAKEN_PROVENANCE="${RUN_DIR}/kraken/${SAMPLE}.Seed${SEED}.mix_generation.txt"
ASM_PROVENANCE="${RUN_DIR}/asm/RRMG_MIX_GENERATION"
ASM_MODE_PROVENANCE="${RUN_DIR}/asm/RRMG_ASSEMBLY_MODE"
ASSEMBLY_MODE=only-assembler
[[ "$ONLY_ASSEMBLER" -eq 1 ]] || ASSEMBLY_MODE=error-correction-and-assembly

require_assembly_mode() {
  marker_matches "$ASM_MODE_PROVENANCE" "$ASSEMBLY_MODE" || {
    echo "[ERROR] Existing assembly mode is different or unrecorded; requested: $ASSEMBLY_MODE" >&2
    echo "        Inspect $ASM_MODE_PROVENANCE and spades.log; rerun Step3 with --force to rebuild." >&2
    return 1
  }
}
echo "[ASSEMBLY] mode=$ASSEMBLY_MODE (fastp QC is controlled separately)"

[[ -n "$TOTAL_PAIRS" ]] || { echo "[ERROR] Plan row has no total_pairs/computed_total_pairs" >&2; exit 1; }
python3 "$FASTQ_MIXER" validate-mix --audit "$MIX_AUDIT" --expected-pairs "$TOTAL_PAIRS" \
  --expected-background-pairs "$BACKGROUND_PAIRS_EXPECTED" --expected-target-pairs "$TARGET_PAIRS_EXPECTED"
MIX_GENERATION="$(awk -F '\t' '$1=="mix_generation_id"{print $2}' "$MIX_AUDIT")"

if [[ "$FORCE" -eq 0 && -s "$SUCCESS_MARKER" && -s "$CONTIGS" && -s "$MAP_SUM" && -s "$MASTER" ]] && \
   awk -F '\t' -v run_dir="$RUN_DIR" 'NR > 1 && $11 == run_dir {found=1} END {exit(found ? 0 : 1)}' "$MASTER"; then
  [[ "$(awk -F '=' '$1=="mix_generation_id"{print $2}' "$SUCCESS_MARKER")" == "$MIX_GENERATION" ]] || {
    echo "[ERROR] Step 3 outputs belong to an older Step 2 read-mix generation: $RUN_DIR" >&2
    echo "        Rerun this row with --force; stale assembly/mapping will not be reused." >&2
    exit 1
  }
  require_assembly_mode
  if [[ "$SKIP_QC" -eq 0 ]]; then
    python3 "$FASTQ_MIXER" validate-fastp --json "$FASTP_JSON" --expected-pairs "$TOTAL_PAIRS"
    marker_matches "$QC_PROVENANCE" "$MIX_GENERATION" || {
      echo "[ERROR] Completed Step 3 has no matching QC provenance marker: $QC_PROVENANCE" >&2
      exit 1
    }
  fi
  if [[ "$SKIP_KRAKEN" -eq 0 ]]; then
    [[ -s "${RUN_DIR}/kraken/${SAMPLE}.Seed${SEED}.report" ]] && \
      marker_matches "$KRAKEN_PROVENANCE" "$MIX_GENERATION" || {
        echo "[ERROR] Completed Step 3 has no matching Kraken2 result for the current mix." >&2
        echo "        Rerun this row with --force to enable/rebuild Kraken2 output." >&2
        exit 1
      }
  fi
  echo "[SKIP] Completed Step 3 result already exists: $RUN_DIR"
  exit 0
fi

# From this point the row is being rebuilt or resumed. A failed attempt must not
# leave an older completion marker that could make the next attempt skip work.
rm -f -- "$SUCCESS_MARKER"

[[ -s "$MIX_R1" && -s "$MIX_R2" ]] || { echo "[ERROR] Mixed reads not found. Run step2 first: $MIX_R1 / $MIX_R2" >&2; exit 1; }
python3 "$FASTQ_MIXER" validate-mix \
  --audit "$MIX_AUDIT" --expected-pairs "$TOTAL_PAIRS" \
  --expected-background-pairs "$BACKGROUND_PAIRS_EXPECTED" --expected-target-pairs "$TARGET_PAIRS_EXPECTED" \
  --r1 "$MIX_R1" --r2 "$MIX_R2"
mkdir -p "$RUN_DIR"/{qc,kraken,asm,map}

QC_R1="$MIX_R1"; QC_R2="$MIX_R2"
if [[ "$SKIP_QC" -eq 0 ]]; then
  QC_R1="${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}_1.fastq.gz"
  QC_R2="${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}_2.fastq.gz"
  QC_REBUILT=0
  if [[ "$FORCE" -eq 0 && -s "$QC_R1" && -s "$QC_R2" ]]; then
    marker_matches "$QC_PROVENANCE" "$MIX_GENERATION" || {
      echo "[ERROR] Existing QC reads belong to an older/unaudited Step 2 mix. Rerun Step 3 with --force." >&2
      exit 1
    }
    echo "[REUSE] Existing QC reads: $QC_R1 / $QC_R2"
  else
    rm -f -- "$QC_R1" "$QC_R2" "$FASTP_JSON" \
      "${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}.fastp.html" "$QC_PROVENANCE"
    run_fastp "$MIX_R1" "$MIX_R2" "$QC_R1" "$QC_R2" "$FASTP_JSON" "${RUN_DIR}/qc/${SAMPLE}.Seed${SEED}.fastp.html"
    QC_REBUILT=1
  fi
  python3 "$FASTQ_MIXER" validate-fastp --json "$FASTP_JSON" --expected-pairs "$TOTAL_PAIRS"
  [[ "$QC_REBUILT" -eq 0 ]] || write_marker "$QC_PROVENANCE" "$MIX_GENERATION"
fi

KRAKEN_RPT="${RUN_DIR}/kraken/${SAMPLE}.Seed${SEED}.report"
KRAKEN_OUT="${RUN_DIR}/kraken/${SAMPLE}.Seed${SEED}.output"
SPECIES_PCT="NA"; SPECIES_CLADE="NA"; SPECIES_TAXON="NA"; SPECIES_TAXID="NA"
if [[ "$SKIP_KRAKEN" -eq 0 ]]; then
  if [[ "$FORCE" -eq 0 && -s "$KRAKEN_RPT" ]]; then
    marker_matches "$KRAKEN_PROVENANCE" "$MIX_GENERATION" || {
      echo "[ERROR] Existing Kraken2 report belongs to an older/unaudited read mix. Rerun Step 3 with --force." >&2
      exit 1
    }
    echo "[REUSE] Existing Kraken2 report: $KRAKEN_RPT"
  else
    rm -f -- "$KRAKEN_RPT" "$KRAKEN_OUT" "$KRAKEN_PROVENANCE"
    run_kraken2 "$QC_R1" "$QC_R2" "$KRAKEN_RPT" "$KRAKEN_OUT"
    write_marker "$KRAKEN_PROVENANCE" "$MIX_GENERATION"
  fi
  if [[ -n "$SPECIES_NAME" ]]; then
    read -r SPECIES_PCT SPECIES_CLADE SPECIES_TAXON SPECIES_TAXID < <(kraken2_species_abundance "$KRAKEN_RPT" "$SPECIES_NAME")
  fi
fi

if [[ "$SKIP_ASSEMBLE" -eq 0 ]]; then
  if [[ "$FORCE" -eq 0 && -s "$CONTIGS" ]]; then
    marker_matches "$ASM_PROVENANCE" "$MIX_GENERATION" || {
      echo "[ERROR] Existing assembly belongs to an older/unaudited read mix. Rerun Step 3 with --force." >&2
      exit 1
    }
    require_assembly_mode
    echo "[REUSE] Existing assembly: $CONTIGS"
  else
    if [[ "$FORCE" -eq 1 ]]; then
      rm -rf -- "${RUN_DIR}/asm"
      mkdir -p "${RUN_DIR}/asm"
    elif [[ -n "$(find "${RUN_DIR}/asm" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
      echo "[ERROR] Incomplete assembly directory exists: ${RUN_DIR}/asm" >&2
      echo "        Inspect it, then rerun this row with --force to replace it." >&2
      exit 1
    fi
    run_metaspades "$QC_R1" "$QC_R2" "${RUN_DIR}/asm"
    write_marker "$ASM_PROVENANCE" "$MIX_GENERATION"
    write_marker "$ASM_MODE_PROVENANCE" "$ASSEMBLY_MODE"
  fi
else
  [[ -s "$CONTIGS" ]] || {
    echo "[ERROR] --skip-assemble requires an existing non-empty file: $CONTIGS" >&2
    exit 1
  }
  marker_matches "$ASM_PROVENANCE" "$MIX_GENERATION" || {
    echo "[ERROR] --skip-assemble cannot reuse contigs from an older/unaudited read mix." >&2
    exit 1
  }
  require_assembly_mode
  echo "[REUSE] Assembly requested with --skip-assemble: $CONTIGS"
fi

echo -e "SampleID\tTargetType\tTargetFA\tContigsFA\tHasHit\tBestAlnLen\tBestIdentity\tTargetBreadthPct\tTargetCoveredBp\tVCF\tVarN" > "$MAP_SUM"
if [[ -n "$GENE_REF" ]]; then contigs_to_target_summary "${SAMPLE}.Seed${SEED}" "$CONTIGS" "$GENE_REF" "GENE" "$MAP_SUM"; fi
if [[ -n "$SPECIES_REF" ]]; then contigs_to_target_summary "${SAMPLE}.Seed${SEED}" "$CONTIGS" "$SPECIES_REF" "SPECIES_GENOME" "$MAP_SUM"; fi
if [[ -z "$GENE_REF" && -z "$SPECIES_REF" ]]; then echo -e "${SAMPLE}.Seed${SEED}\tNA\tNA\t${CONTIGS}\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "$MAP_SUM"; fi

mkdir -p result

# PBS array tasks share one summary file. Serialize updates and replace an
# earlier row for the same run directory instead of accumulating duplicates.
(
  command -v flock >/dev/null 2>&1 || {
    echo "[ERROR] flock is required for safe parallel summary updates" >&2
    exit 1
  }
  flock -x 9
  MASTER_TMP="${MASTER}.tmp.${PBS_JOBID:-manual}.${ROW_INDEX}.$$"
  if [[ -s "$MASTER" ]]; then
    awk -F '\t' -v run_dir="$RUN_DIR" 'NR == 1 || $11 != run_dir' "$MASTER" > "$MASTER_TMP"
  else
    echo -e "Seed\tSample\tTargetAbundance\tTotalPairs\tKrakenSpecies\tSpeciesPct\tSpeciesCladeReads\tSpeciesTaxonReads\tSpeciesTaxID\tMapSummary\tRunDir" > "$MASTER_TMP"
  fi
  echo -e "${SEED}\t${SAMPLE}\t${TARGET_ABUNDANCE}\t${TOTAL_PAIRS}\t${SPECIES_NAME:-NA}\t${SPECIES_PCT}\t${SPECIES_CLADE}\t${SPECIES_TAXON}\t${SPECIES_TAXID}\t${MAP_SUM}\t${RUN_DIR}" >> "$MASTER_TMP"
  mv -f -- "$MASTER_TMP" "$MASTER"
) 9>"${MASTER}.lock"

printf 'run_id=%s\nrow_index=%s\nmix_generation_id=%s\nassembly_mode=%s\ncompleted=%s\n' \
  "$RUN_ID" "$ROW_INDEX" "$MIX_GENERATION" "$ASSEMBLY_MODE" "$(date -Is)" > "$SUCCESS_MARKER"
echo "[DONE] Downstream summary -> $MASTER"
