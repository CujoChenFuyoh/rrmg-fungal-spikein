#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  rrmg_step2_make_target_and_merge.sh plan.tsv ROW_INDEX -w /path/to/project \
    --background-index temp/rrmg_background/background_index.tsv

Purpose:
  For one plan row / seed, generate ONLY target-species reads with that seed,
  then merge them with precomputed background reads.

Inputs:
  permutation plan, background_index.tsv, background FASTQ, target genome FASTA

Outputs:
  <run_dir>/reads/<sample>.Seed<seed>_{1,2}.fastq.gz
  <run_dir>/design/merge_design.tsv

Important:
  Use --keep-target-reads when Step 5 or the target-only CYP51A comparison will
  be run later. Without it, target-only FASTQs are deleted after merging.

Required environment:
  RANDOMREADSMG_SH      Absolute path to randomreadsmg.sh
USAGE
}

PLAN="${1:-}"
ROW_INDEX="${2:-${PBS_ARRAY_INDEX:-}}"
if [[ -z "$PLAN" || -z "$ROW_INDEX" ]]; then usage; exit 1; fi
shift 2 || true

WD=""
BACKGROUND_INDEX="temp/rrmg_background/background_index.tsv"
FORCE=0
KEEP_TARGET=0
READLEN_DEFAULT=150

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WD="$2"; shift 2;;
    --background-index) BACKGROUND_INDEX="$2"; shift 2;;
    --force) FORCE=1; shift 1;;
    --keep-target-reads) KEEP_TARGET=1; shift 1;;
    --readlen) READLEN_DEFAULT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

: "${WD:?Missing -w}"
: "${RANDOMREADSMG_SH:?RANDOMREADSMG_SH is not set}"
WD="$(mkdir -p "$WD" && cd "$WD" && pwd)"
cd "$WD"
[[ -s "$PLAN" ]] || { echo "[ERROR] Plan not found: $PLAN" >&2; exit 1; }
[[ -s "$BACKGROUND_INDEX" ]] || { echo "[ERROR] Background index not found: $BACKGROUND_INDEX" >&2; exit 1; }
[[ -s "$RANDOMREADSMG_SH" ]] || { echo "[ERROR] RANDOMREADSMG_SH not found: $RANDOMREADSMG_SH" >&2; exit 1; }

ROW_ENV="$(mktemp "${WD}/temp/rrmg_row_env.XXXXXX")"
python3 - "$PLAN" "$ROW_INDEX" "$WD" "$BACKGROUND_INDEX" "$READLEN_DEFAULT" > "$ROW_ENV" <<'PY'
import csv, gzip, hashlib, math, re, shlex, sys
from pathlib import Path

plan = Path(sys.argv[1]); row_index = int(sys.argv[2]); wd = Path(sys.argv[3]); bg_index = Path(sys.argv[4]); readlen_default = int(float(sys.argv[5]))

def safe(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip())
    return s.strip('_') or 'NA'

def norm(s): return re.sub(r"[\s_\-]+", "", str(s)).lower()
def resolve(p):
    if not p: return ''
    q = Path(p)
    return str(q if q.is_absolute() else wd/q)
def q(v): return shlex.quote(str(v))
def open_text(path):
    path = Path(path)
    return gzip.open(path,'rt',errors='ignore') if str(path).endswith('.gz') else path.open('r',errors='ignore')
def fa_len(path):
    n=0
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith('>'): n += len(line.strip())
    return n

def read_manifest(manifest):
    with Path(manifest).open(newline='') as fh:
        sample=fh.read(8192); fh.seek(0)
        try: dialect=csv.Sniffer().sniff(sample, delimiters='\t,;')
        except csv.Error: dialect=csv.excel_tab
        rows=list(csv.DictReader(fh,dialect=dialect))
    out=[]
    for i,r in enumerate(rows,1):
        keys={norm(k):k for k in r}
        def get(*names):
            for name in names:
                if norm(name) in keys: return r.get(keys[norm(name)],'')
            return ''
        out.append({'gid':get('genome_ID','genome','id') or f'Genome{i}', 'name':get('name','species','taxon'), 'taxid':get('taxid','NCBI ID','NCBI_ID','TaxID'), 'fasta':resolve(get('fasta','path_to_fasta_file','genome_path','path')), 'role':(get('role') or 'background').lower()})
    return out

def is_target(r, target_name, target_taxid, target_gid):
    return (r.get('role')=='target' or (target_gid and r.get('gid')==target_gid) or (target_taxid and target_taxid not in {'NA','0',''} and str(r.get('taxid'))==str(target_taxid)) or (target_name and r.get('name','').strip().lower()==target_name.strip().lower()))

with plan.open(newline='') as fh:
    reader=csv.DictReader(fh, delimiter='\t')
    chosen=None
    for row in reader:
        if int(row.get('run_index','0') or 0)==row_index:
            chosen=row; break
if chosen is None: raise SystemExit(f'[ERROR] Row index not found: {row_index}')

manifest=resolve(chosen.get('manifest',''))
ta=float(chosen.get('target_abundance') or 0)
readlen=int(float(chosen.get('readlen') or readlen_default))
total_pairs=int(round(float(chosen.get('total_pairs') or 0)))
if total_pairs <= 0:
    total_bases=chosen.get('total_bases','')
    total_cov=chosen.get('total_cov','')
    csize=chosen.get('community_size_bp','')
    if total_bases:
        total_pairs=int(round(float(total_bases)/(readlen*2)))
    elif total_cov and csize:
        total_pairs=int(round(float(total_cov)*float(csize)/(readlen*2)))
    else:
        raise SystemExit('[ERROR] Step2 currently requires total_pairs, total_bases, or total_cov+community_size_bp in plan')
target_pairs=int(round(total_pairs*ta))
bg_pairs=total_pairs-target_pairs
sample=chosen.get('sample') or chosen.get('run_id') or f'row{row_index}'
run_id=chosen.get('run_id') or f'{sample}__Seed{chosen.get("seed",row_index)}'
seed=int(float(chosen.get('seed') or row_index))
target_seed=seed + 2900000

target_name=chosen.get('target_name','')
target_taxid=chosen.get('target_taxid','')
target_gid=chosen.get('target_genome_id','')
target_fasta=resolve(chosen.get('target_fasta',''))
if not target_fasta:
    for r in read_manifest(manifest):
        if is_target(r,target_name,target_taxid,target_gid):
            target_fasta=r['fasta']; break
if target_pairs > 0 and (not target_fasta or not Path(target_fasta).exists()):
    raise SystemExit(f'[ERROR] Target FASTA not found: {target_fasta}')
target_len=fa_len(target_fasta) if target_fasta and Path(target_fasta).exists() else 0
target_cov=(target_pairs*readlen*2/target_len) if target_len else 0

bg_key_base=safe(sample)
digest=hashlib.md5(f"{manifest}|{ta}|{total_pairs}|{readlen}|{chosen.get('abundance_mode','read_fraction') or 'read_fraction'}".encode()).hexdigest()[:8]
bg_key=f"{bg_key_base}__bg_{digest}"
with bg_index.open(newline='') as fh:
    rows=list(csv.DictReader(fh, delimiter='\t'))

def _as_float(x, default=None):
    try:
        if x is None or str(x).strip() == '':
            return default
        return float(str(x).strip())
    except Exception:
        return default

def _as_int(x, default=None):
    y = _as_float(x, None)
    if y is None:
        return default
    return int(round(y))

def _norm_path(x):
    if x is None:
        return ''
    s = str(x).strip()
    if not s:
        return ''
    # Resolve relative paths against WD, but tolerate Imperial /user/ vs /ephemeral/user/ aliases.
    q = Path(s)
    if not q.is_absolute():
        q = wd / q
    s = str(q)
    s = s.replace('/rds/general/ephemeral/user/', '/rds/general/user/')
    return s.rstrip('/')

def _float_same(a, b, tol=1e-9):
    aa = _as_float(a, None)
    bb = _as_float(b, None)
    return (aa is not None) and (bb is not None) and abs(aa - bb) <= tol

match=[r for r in rows if r.get('bg_key')==bg_key]

# Robust fallback: older/newer scripts may create different bg_key strings because of
# sample-name formatting such as pairs100000000 vs pairs100000000.0, or path aliases.
# Match the actual design fields instead of the digest when exact key matching fails.
if not match:
    manifest_norm = _norm_path(manifest)
    candidates = []
    for r in rows:
        r_manifest = _norm_path(r.get('manifest', ''))
        r_ta = r.get('target_abundance', '')
        r_pairs = r.get('total_pairs', '')
        r_readlen = r.get('readlen', '')
        if r_manifest != manifest_norm:
            continue
        if not _float_same(r_ta, ta, tol=1e-8):
            continue
        if _as_int(r_pairs, -1) != int(total_pairs):
            continue
        if _as_int(r_readlen, -1) != int(readlen):
            continue
        candidates.append(r)
    if len(candidates) == 1:
        match = candidates
        print(f'[WARN] Exact background key not found, but matched by design fields: {candidates[0].get("bg_key", "")}', file=sys.stderr)
    elif len(candidates) > 1:
        keys = ','.join([c.get('bg_key','') for c in candidates[:5]])
        raise SystemExit(f'[ERROR] Multiple background rows match design fields for {bg_key}: {keys}')

if not match:
    nearby = []
    for r in rows:
        if _float_same(r.get('target_abundance',''), ta, tol=1e-8) and _as_int(r.get('total_pairs',''), -1) == int(total_pairs):
            nearby.append(r.get('bg_key',''))
    msg = [
        f'[ERROR] Background key not found in index: {bg_key}',
        f'        Wanted manifest={manifest}',
        f'        Wanted target_abundance={ta}',
        f'        Wanted total_pairs={total_pairs}',
        f'        Wanted readlen={readlen}',
        f'        Background index rows={len(rows)}',
    ]
    if nearby:
        msg.append('        Nearby bg_key with same abundance/pairs: ' + ', '.join(nearby[:5]))
    else:
        msg.append('        No background row found with same target_abundance and total_pairs.')
    raise SystemExit('\n'.join(msg))
bg=match[0]
bg_r1=bg['bg_r1']; bg_r2=bg['bg_r2']
if not Path(bg_r1).exists() or not Path(bg_r2).exists():
    raise SystemExit(f'[ERROR] Background FASTQ missing: {bg_r1} / {bg_r2}')

run_dir=wd/'temp'/'rrmg_permutation'/safe(run_id)
reads_dir=run_dir/'reads'; design_dir=run_dir/'design'
reads_dir.mkdir(parents=True, exist_ok=True); design_dir.mkdir(parents=True, exist_ok=True)
run_tag=f"Seed{seed}"
target_r1=reads_dir/f"{safe(target_gid or 'target')}.{run_tag}.target_1.fastq.gz"
target_r2=reads_dir/f"{safe(target_gid or 'target')}.{run_tag}.target_2.fastq.gz"
mix_r1=reads_dir/f"{safe(sample)}.{run_tag}_1.fastq.gz"
mix_r2=reads_dir/f"{safe(sample)}.{run_tag}_2.fastq.gz"
meta=design_dir/'merge_design.tsv'
with meta.open('w') as fh:
    fh.write('key\tvalue\n')
    for k,v in [('row_index',row_index),('run_id',run_id),('sample',sample),('bg_key',bg_key),('background_r1',bg_r1),('background_r2',bg_r2),('target_name',target_name),('target_taxid',target_taxid),('target_genome_id',target_gid),('target_fasta',target_fasta),('target_abundance',ta),('total_pairs',total_pairs),('background_pairs',bg_pairs),('target_pairs',target_pairs),('target_seed',target_seed),('target_genome_bp',target_len),('target_cov',target_cov),('mixed_r1',mix_r1),('mixed_r2',mix_r2),('readlen',readlen)]:
        fh.write(f'{k}\t{v}\n')

for k,v in {
'RUN_ID':run_id,'SAMPLE':sample,'READLEN':readlen,'TARGET_FASTA':target_fasta,'TARGET_PAIRS':target_pairs,'TARGET_SEED':target_seed,'TARGET_R1':target_r1,'TARGET_R2':target_r2,'BG_R1':bg_r1,'BG_R2':bg_r2,'MIX_R1':mix_r1,'MIX_R2':mix_r2,'RUN_DIR':run_dir,'TARGET_COV':target_cov,'TARGET_GENOME_BP':target_len}.items():
    print(f'{k}={q(v)}')
PY
# shellcheck disable=SC1090
source "$ROW_ENV"
rm -f "$ROW_ENV"

mkdir -p "$(dirname "$TARGET_R1")" "$(dirname "$MIX_R1")"
if [[ "$FORCE" -eq 0 && -s "$MIX_R1" && -s "$MIX_R2" ]]; then
  echo "[MERGE] exists, skip: $MIX_R1 / $MIX_R2"
  exit 0
fi

if [[ "$TARGET_PAIRS" -gt 0 ]]; then
  echo "[TARGET-SIM] sample=$SAMPLE target_pairs=$TARGET_PAIRS target_cov=$TARGET_COV seed=$TARGET_SEED"
  bash "$RANDOMREADSMG_SH" \
    -Xmx"${BBMAP_XMX:-4g}" t="${THREADS:-1}" overwrite=t \
    ref="$TARGET_FASTA" out="$TARGET_R1" out2="$TARGET_R2" \
    reads="$TARGET_PAIRS" paired=t minlength="$READLEN" maxlength="$READLEN" seed="$TARGET_SEED" \
    2> "${TARGET_R1}.randomreadsmg.log"
  [[ -s "$TARGET_R1" && -s "$TARGET_R2" ]] || { echo "[ERROR] target randomreadsmg failed; see ${TARGET_R1}.randomreadsmg.log" >&2; exit 1; }
  rm -f "$MIX_R1" "$MIX_R2"
  cat "$BG_R1" "$TARGET_R1" > "$MIX_R1"
  cat "$BG_R2" "$TARGET_R2" > "$MIX_R2"
  if [[ "$KEEP_TARGET" -eq 0 ]]; then
    rm -f "$TARGET_R1" "$TARGET_R2"
  fi
else
  echo "[TARGET-SIM] target_pairs=0; mixed reads are symlinks to background reads"
  rm -f "$MIX_R1" "$MIX_R2"
  ln -s "$BG_R1" "$MIX_R1"
  ln -s "$BG_R2" "$MIX_R2"
fi

[[ -s "$MIX_R1" && -s "$MIX_R2" ]] || { echo "[ERROR] merge failed: $MIX_R1 / $MIX_R2" >&2; exit 1; }
echo "[DONE] Mixed reads -> $MIX_R1 / $MIX_R2"
