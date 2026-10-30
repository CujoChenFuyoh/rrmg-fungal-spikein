#!/usr/bin/env bash
# RRMG 0.3.4 本地 Ubuntu/WSL 逐段运行手册
# ------------------------------------------------------------
# 用法：在编辑器打开本文件，按编号依次复制到同一个 Ubuntu Bash 终端。
# 带反斜杠的多行命令、if/fi、Python heredoc 必须作为一个完整块执行。
# 不在 Windows PowerShell 里执行，不需要 qsub，不调用 PBS 文件。
# 建议首次逐段检查，不直接 bash runbook.sh 一次跑到底。
# 在仓库目录开始；先设置下方输入路径。安装细节见 docs/INSTALL.zh-CN.md。
#
# 测试设计：2 个演示背景物种 + Af293；1 个丰度/深度条件；Seed 1、2。
# 设计：每重复 1,000,000 对 reads，目标 reads 比例 0.01，read length 150。
# 这是执行检查，低深度未检出 CYP51A 也可能是正常结果。
# 当前 Step2 强制 target/background/total reads 与设计完全一致，任一不符即停止。
# 不默认运行 Kraken2/VAMB；不做突变设计、变异判定或自动删除 reads。
#
# 失败时停在当前步骤，先看报错/日志，不继续粘贴下一步骤。
# 续跑：执行第 0 段恢复变量/环境，再 source "$RRMG_CONFIG"；
#       跳过 init、manifest 和 plan 生成，直接运行尚未完成的 stage。
# 不在已有服务器实验目录中执行第 1–5 段。


# ===== 0. 路径和环境：先核对下面的路径 =====
CODE="${CODE:-$PWD}"
export WD="${WD:-$HOME/rrmg_local_test_01}"
BBMAP_DIR_LOCAL="${BBMAP_DIR_LOCAL:-}"
# 环境内BBMap在激活后自动解析；外部安装时先 export BBMAP_DIR_LOCAL=/absolute/path/bbmap。
GENE_SOURCE="${GENE_SOURCE:-}"
# 先 export GENE_SOURCE=/absolute/path/CYP51A_wt.fasta；必须是已核验的基因参考。
# runbook、R 绘图脚本和辅助脚本均随仓库提供，统一通过上面的 CODE 调用。
SITES_SOURCE="${SITES_SOURCE:-}"
CATALOGUE_SOURCE="${CATALOGUE_SOURCE:-}"
# 位点注释可选；需要时同时设置以上两个TSV。Windows R可设置WIN_RSCRIPT/WIN_R_LIB。
export RRMG_CONFIG="$WD/rrmg_config.sh"

CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
source "$CONDA_SH"
conda activate fungamr_meta
if [[ -z "$BBMAP_DIR_LOCAL" ]] && command -v randomreadsmg.sh >/dev/null; then
    BBMAP_DIR_LOCAL="$(dirname "$(readlink -f "$(command -v randomreadsmg.sh)")")"
fi

# 报错即停。不要把报错后的命令也粘贴进去。
set -eo pipefail
cd "$CODE"
test -s VERSION
: "${GENE_SOURCE:?Set GENE_SOURCE to the verified target-gene FASTA}"
: "${BBMAP_DIR_LOCAL:?Set BBMAP_DIR_LOCAL or install BBMap in the active environment}"
test -s "$BBMAP_DIR_LOCAL/randomreadsmg.sh"
test -s "$BBMAP_DIR_LOCAL/randomreads.sh"
test -s "$GENE_SOURCE"

# datasets 和 Java 都来自已激活的 fungamr_meta，不切换下载环境。
command -v python
command -v datasets
datasets version
command -v java
readlink -f "$(command -v java)"
java -version
df -h "$HOME"
free -h


# ===== 1. 安装本地最新代码，创建独立测试配置（只做一次） =====
# 从当前检出的源码安装。
python -m pip install .
rrmg --version

# 防止误把新实验写进已有目录；已存在时换一个 WD 或按开头说明续跑。
if [[ -e "$WD" ]]; then
    printf '[STOP] 测试目录已存在，请换 WD 或跳过初始化：%s\n' "$WD" >&2
    exit 1
fi

rrmg init \
    --workdir "$WD" \
    --pipeline-dir "$CODE" \
    --bbmap-dir "$BBMAP_DIR_LOCAL"

source "$RRMG_CONFIG"
mkdir -p "$WD/input/rrmg" "$WD/ref" "$WD/temp" "$WD/result"

# 只复制现有基因参考到新测试目录，不改变原文件，不自动重新提取基因。
cp -- "$GENE_SOURCE" "$WD/ref/CYP51A_wt.fasta"
cmp -- "$GENE_SOURCE" "$WD/ref/CYP51A_wt.fasta"
head -n 2 "$WD/ref/CYP51A_wt.fasta"

# 可选：nano "$RRMG_CONFIG"
# 修改配置后必须重新执行：source "$RRMG_CONFIG"
# 当前默认：Step1 4线程/4g heap，Step2 1线程/4g heap；
# Step3 4线程/16GB SPAdes；Step4/5 4线程；LOCAL_JOBS=1。
# Step3 默认 STEP3_ONLY_ASSEMBLER=1：保留 fastp，但跳过 BayesHammer reads 纠错。
# 先检查机器剩余内存。资源不足时在配置中调低线程/内存或减少测试深度。
# 下方提交命令用 -t/--threads 临时设置每行线程；Step1–3 还可用
# -m/--memory 临时设置每行程序内存（正整数 GB，不设人为上限），无需改配置。
# -m 不会向 PBS 申请资源；PBS mem/walltime 仍必须在用户自己的作业脚本中设置。
# 保持 STEP2_KEEP_TARGET_READS=1、RUN_KRAKEN=0；配置中的各 STEP*_FORCE 保持 0。
# 下方 Step1–5 示例显式加入 --force，供独立测试目录从头重建结果；
# 恢复或续跑正式数据时应删除 --force，让完整的已有输出被安全复用。
# Java 默认留空，使用激活环境的 PATH，不需要写死 JAVA_HOME。
rrmg doctor --config "$RRMG_CONFIG"
# 此时 plan/background/genome 尚未生成，相关 WARN 正常；工具 ERROR 必须处理。


# ===== 2. Prepare manifest：从演示物种表生成下载清单 =====
# 演示表包含两个背景物种，不是原研究的群落设计。
column -t -s $'\t' "$CODE/examples/taxa_demo.tsv"
test ! -e "$WD/input/rrmg/community_manifest.tsv"

python3 "$CODE/rrmg_prepare_manifest.py" \
    --input "$CODE/examples/taxa_demo.tsv" \
    --outdir "$WD/input/rrmg" \
    --topN 2 \
    --target-name 'Aspergillus fumigatus' \
    --target-taxid 330879 \
    --target-genome-id Afum \
    --target-fasta "$WD/ref/Afum_Af293.fasta"

column -t -s $'\t' "$WD/input/rrmg/community_manifest.tsv"
# 预期：3 条数据行；2 background + 1 target；target 的 rel_abund 为 0。
# Af293 下载使用菌株 taxid；目标比例在后续 plan 中指定。


# ===== 3. 下载 genome：第一次需要网络，可能耗时 =====
# 保留原始下载报告以核对 assembly 身份；已存在的非空 FASTA 会被复用。
# 下载器按物种查 reference，并可能回退到较大的 assembly 数据包。
# --threads 目前不控制下载并行，因此这里不设置。
bash "$CODE/rrmg_download_genomes.sh" \
    --manifest "$WD/input/rrmg/community_manifest.tsv" \
    --out-manifest "$WD/input/rrmg/community_manifest.downloaded.tsv" \
    --keep-tmp

cat "$WD/input/rrmg/download_failed.tsv"
column -t -s $'\t' "$WD/input/rrmg/community_manifest.downloaded.tsv"
test -s "$WD/ref/Afum_Af293.fasta"

# 此演示要求全部三个 genome 下载成功，避免悄悄接受变化后的背景。
# 如果这里失败，先检查 input/rrmg/download_logs/，不要继续生成 plan。
test ! -s "$WD/input/rrmg/download_failed.tsv"
test "$(awk 'NR>1 {n++} END {print n+0}' "$WD/input/rrmg/community_manifest.downloaded.tsv")" -eq 3


# ===== 4. 用少量模拟 reads 检查 genome/BBMap/Java 兼容性 =====
source "$CODE/rrmg_environment.sh"
export THREADS=1
export BBMAP_XMX=2g

bash "$CODE/rrmg_test_randomreadsmg_genomes.sh" \
    --manifest "$WD/input/rrmg/community_manifest.downloaded.tsv" \
    --out-manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
    -w "$WD" \
    --test-pairs 200 --readlen 150 --seed 999

column -t -s $'\t' "$WD/input/rrmg/randomreadsmg_genome_test/randomreadsmg_test_summary.tsv"
test ! -s "$WD/input/rrmg/randomreadsmg_genome_test/randomreadsmg_test_failed.tsv"
test "$(awk 'NR>1 {n++} END {print n+0}' "$WD/input/rrmg/community_manifest.tested.tsv")" -eq 3

# 固定这次测试的输入和环境记录；不会重新创建 Conda 环境。
sha256sum "$WD/ref/"*.fasta "$WD/input/rrmg/genomes/"*.fna > "$WD/result/reference_files.sha256"
conda list --explicit > "$WD/result/conda_explicit.txt"
cp -- "$CODE/VERSION" "$WD/result/workflow_VERSION.txt"
git -C "$CODE" rev-parse HEAD > "$WD/result/workflow_commit.txt"
git -C "$CODE" status --short > "$WD/result/workflow_status.txt"
# 本地未提交修改不会由 commit SHA 代表，因此另存实际核心脚本快照。
rrmg export "$WD/code_snapshot_0.3.4"


# ===== 5. 生成测试 plan 和背景设计（只做一次） =====
test ! -e "$WD/input/rrmg/rrmg_permutation_plan.tsv"
test ! -e "$WD/temp/rrmg_background/background_plan.tsv"

python3 "$CODE/rrmg_make_plan.py" \
    --manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
    --out "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
    --sample-prefix LocalPilot_Afum \
    --target-name 'Aspergillus fumigatus' \
    --target-taxid 330879 \
    --target-genome-id Afum \
    --target-fasta "$WD/ref/Afum_Af293.fasta" \
    --target-abundances '0.01' \
    --total-pairs '1000000' \
    --seeds '1-2' \
    --readlen 150

python3 "$CODE/rrmg_make_background_plan.py" \
    --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
    -w "$WD" \
    --outdir "$WD/temp/rrmg_background" \
    --background-seed 1000

wc -l "$WD/input/rrmg/rrmg_permutation_plan.tsv" "$WD/temp/rrmg_background/background_plan.tsv"
# 预期分别 3 行、2 行（都包含表头），即两个 permutation、一个背景条件。
rrmg doctor --config "$RRMG_CONFIG"


# ===== 6. Step1：只生成一次共享背景 reads =====
# 设置参数（直接修改下面命令，无需编辑 config）：
#   -t 4 / --threads 4      # 每个 background 条件使用 4 线程
#   -m 4 / --memory 4       # 每个条件给 BBMap Java 4 GB heap；只检查为正整数
# 仍由 config 控制：BACKGROUND_SEED=1000；修改 seed 后必须重跑依赖它的 Step2–5。
# 省略 -t/-m 时才使用 STEP1_THREADS/STEP1_BBMAP_XMX 的配置值。
# 提交范围：step1 后的 1 是 background_plan.tsv 中的 bg_index。
# 并发设置：--jobs 1 表示一次运行 1 个背景条件；总 CPU/内存近似按 jobs 倍增。
# PBS 包装脚本申请的 ncpus 至少应覆盖 --jobs × 这里的 -t。
# 数字 1 是 background_plan 的 bg_index，不是 permutation 的 run_index。
bash "$CODE/rrmg_submit.sh" step1 1 -t 4 -m 4 --dry-run
bash "$CODE/rrmg_submit.sh" step1 1 -t 4 -m 4 --jobs 1
find "$WD/temp/rrmg_background" -name '*.pair_count.tsv' -type f -size +0c -print \
    -exec column -t -s $'\t' {} \;
# source_pairs 可大于设计背景数；Step2 会均匀下采样到精确 background_pairs。


# ===== 7. Step2：两个重复分别模拟 target reads，再与背景合并 =====
# 设置参数（直接修改下面命令，无需编辑 config）：
#   -t 1 / --threads 1      # 建议保持 1，使固定 seed 的 BBMap 模拟可重复
#   -m 4 / --memory 4       # 每个 run_index 给 BBMap Java 4 GB heap
# 仍由 config 控制：STEP2_KEEP_TARGET_READS=1；Step5/alignment-block 需要它。
# 省略 -t/-m 时才使用 STEP2_THREADS/STEP2_BBMAP_XMX 的配置值。
# 提交范围：1-2 来自 rrmg_permutation_plan.tsv 的 run_index，可改为 1、1-20 等。
# 并发设置：--jobs 1 表示一次运行 1 个样品；提高后按 jobs 倍增 CPU、heap 和 I/O。
# 强制重做 Step2 会产生新的 read generation；同一批样品的 Step3–5 也必须重跑。
# 这里 1-2 是 permutation plan 的 run_index，不是文本行号或物种 ID。
bash "$CODE/rrmg_submit.sh" step2 1-2 -t 1 -m 4 --dry-run
bash "$CODE/rrmg_submit.sh" step2 1-2 -t 1 -m 4 --jobs 1
find "$WD/temp/rrmg_permutation" -path '*/design/read_mix_audit.tsv' -type f -size +0c -print \
    -exec column -t -s $'\t' {} \;
# 每个 audit 必须是 gzip_members_per_output=1，且 background+target=total。


# ===== 8. Step3：fastp QC → metaSPAdes → contig-to-gene mapping =====
# 设置参数（直接修改下面命令，无需编辑 config）：
#   -t 4 / --threads 4      # 每个样品供 fastp/metaSPAdes 使用的线程数
#   -m 16 / --memory 16     # 每个样品传给 metaSPAdes 的内存值（GB，无人为上限）
# 仍由 config 控制：
#   STEP3_PHRED_OFFSET=33       # 本流程模拟 Illumina reads 的质量编码
#   STEP3_ONLY_ASSEMBLER=1      # 默认：fastp 后以 --only-assembler 拼接，跳过 BayesHammer
#                              # 改为 0 则恢复原生 reads 纠错；不要混用两种模式做同组比较
#   STEP3_SKIP_QC=0             # 设1使用raw reads跳过fastp；不是复用QC的开关
#   STEP3_SKIP_ASSEMBLY=0       # 设 1 仅在已有可验证 assembly 时跳过拼接
#   RUN_KRAKEN=0                # 设 1 时还必须设置 KRAKEN_DB
# 省略 -t/-m 时才使用 STEP3_THREADS/STEP3_SPADES_MEM_GB 的配置值。
# 并发设置：--jobs 1 表示一次拼接 1 个样品；PBS ncpus 至少为 jobs × -t，
#           PBS mem 应覆盖 jobs × -m，并额外预留其他进程和系统开销。
# 这个 QC 是 fastp，不是 FastQC。不调用 Kraken2，不进行 VAMB binning。
# 老 config 没有 STEP3_ONLY_ASSEMBLER 时也默认 1；需要纠错时在 RRMG_CONFIG 文件中设为 0。
# 修改模式后需 --force 重建 Step3；仅 source 后设置未写入 config 的变量不会覆盖配置文件。
# 日志应显示 [ASSEMBLY] mode=only-assembler；原生命令应包含 --only-assembler。
# 完成后 asm/RRMG_ASSEMBLY_MODE 和 STEP3_SUCCESS 记录本次模式。
# 已有旧 assembly 若缺少模式记录会要求重建，不会将旧 contigs 静默当成新模式结果。
# 下面 -t/-m/--jobs 是用户可修改的示例资源值，不是流程对服务器资源的限制。
bash "$CODE/rrmg_submit.sh" step3 1-2 -t 4 -m 16 --jobs 1
find "$WD/temp/rrmg_permutation" -path '*/asm/contigs.fasta' -type f -size +0c -print
# 应看到两个 contigs.fasta。无目标基因命中不等于 assembly 运行失败。
# 重新拼接后继续第 11–14 段，重收集汇总并重建 figure_data；不能只重画旧汇总。


# ===== 9. Step4：mixed QC + target-only reads → CYP51A =====
# 设置参数（线程直接写在下面命令；无需编辑 config）：
#   -t 4 / --threads 4          # 每个样品的 minimap2/samtools 线程数
#   GENE_REF="ref/CYP51A_wt.fasta"  # 相对 WD 或绝对路径；必须与后续绘图参考一致
# Step4 没有可靠的程序总内存参数，因此不接受 -m；PBS mem 在作业脚本中按服务器设置。
# 提交范围和 --jobs 规则与 Step2 相同；PBS ncpus 至少为 jobs × -t。
# Step4 同时使用 Step3 QC reads 和 Step2 保留的 target-only reads；分别生成
# BAM、逐碱基 depth、coverage summary 和 read-level PAF。旧 mixed-read 结果完整时，
# 新版只补 target-only mapping；--force 才会重建两类 Step4 输出。
bash "$CODE/rrmg_submit.sh" step4 1-2 -t 4 --jobs 1
# 默认保留PAF用于read-block图；不使用 --no-paf。
# Step3 自动增加 contig cs/cg PAF；Step4 自动生成两类 reads 的 raw/Q20 碱基计数。
# 已有旧项目只补证据（不重做 reads/QC/assembly/read mapping）：
# bash "$CODE/rrmg_backfill_cyp51a_evidence.sh" -w "$WD" -r all -t 4


# ===== 10. Step5：target-only reads → Af293 genome coverage =====
# 设置参数（线程直接写在下面命令；无需编辑 config）：
#   -t 4 / --threads 4                # 每个样品的 minimap2/samtools 线程数
#   TARGET_GENOME="ref/Afum_Af293.fasta"  # 相对 WD 或绝对路径
#   STEP45_MAKE_BINS=1                   # 生成固定窗口 coverage 表
#   STEP45_BIN_SIZE=10000                # coverage 窗口大小（bp），不是 MAG bin
# Step5 没有可靠的程序总内存参数，因此不接受 -m；PBS mem 在作业脚本中按服务器设置。
# Step5 必须有 STEP2_KEEP_TARGET_READS=1 生成的 target-only 双端 reads。
# 提交范围和 --jobs 规则与 Step2 相同；提高 jobs 时 CPU、内存和磁盘 I/O 同时增加。
bash "$CODE/rrmg_submit.sh" step5 1-2 -t 4 --jobs 1
# step45 1-2 可代替上面两步，但本手册拆开以便定位失败位置。
# Step5 的 coverage windows 不是 genome binning 结果。


# ===== 11. 汇总绘图输入到 result/ =====
python3 "$CODE/rrmg_collect_map_results.py" \
    --summary "$WD/result/rrmg_split_pipeline_summary.tsv" \
    --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
    --root "$WD/temp/rrmg_permutation" \
    --out "$WD/result/rrmg_target_gene_map_merged.tsv" \
    --wd "$WD"

bash "$CODE/rrmg_collect_supplementary_downstream.sh" \
    "$WD" "$WD/temp/rrmg_permutation"


# ===== 12. 核对样品 ID、必需输出和 target reads 实际数量 =====
# 下面整个 Python 块一次复制，包括最后的 PY；只读检查，不修改结果。
# 预期每张汇总表恰好包含 plan 的两个 RunID；缺失不能当作生物学零值。
python3 - "$WD" <<'PY'
import csv
import json
import sys
from pathlib import Path

wd = Path(sys.argv[1])
def rows(path):
    with path.open(newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

plan = rows(wd / 'input/rrmg/rrmg_permutation_plan.tsv')
expected = {r['run_id'] for r in plan}
assert len(expected) == len(plan), 'Duplicate run_id in plan'
print('Expected samples:', len(expected))
required = ['design/read_mix_audit.tsv', 'asm/contigs.fasta', 'map/target_call_summary.tsv',
            'readmap/reads_to_GENE.coverage_summary.tsv', 'readmap/reads_to_GENE.depth.tsv',
            'readmap/reads_to_GENE.paf.gz',
            'readmap/target_only_reads_to_GENE.coverage_summary.tsv',
            'readmap/target_only_reads_to_GENE.depth.tsv',
            'readmap/target_only_reads_to_GENE.paf.gz',
            'readmap/reads_to_GENE.base_counts.raw.tsv.gz',
            'readmap/reads_to_GENE.base_counts.q20.tsv.gz',
            'readmap/target_only_reads_to_GENE.base_counts.raw.tsv.gz',
            'readmap/target_only_reads_to_GENE.base_counts.q20.tsv.gz',
            'target_readmap/afum_target_reads_to_genome.summary.tsv',
            'target_readmap/afum_target_reads_to_genome.bins.tsv']
missing_files = []
for relative in required:
    found = 0
    for run_id in sorted(expected):
        path = wd / 'temp/rrmg_permutation' / run_id / relative
        if path.is_file() and path.stat().st_size > 0:
            found += 1
        else:
            missing_files.append(str(path))
    print(f'{relative}: {found}/{len(expected)}')
for name in ['rrmg_target_gene_map_merged.tsv', 'rrmg_reads_to_gene_depth_summary.tsv',
             'rrmg_target_only_reads_to_gene_depth_summary.tsv',
             'rrmg_afum_target_reads_genome_coverage_summary.tsv']:
    table = rows(wd / 'result' / name)
    ids = [r['RunID'] for r in table]
    print(f'{name}: {len(ids)} rows; missing={sorted(expected-set(ids))}; extra={sorted(set(ids)-expected)}')
    assert set(ids) == expected and len(ids) == len(expected), 'Incomplete/duplicate summary: ' + name
assert not missing_files, 'Missing output files:\n' + '\n'.join(missing_files)
# cs PAF可为空（没有命中），但完成侧录必须非空。
for run_id in sorted(expected):
    cs = wd / 'temp/rrmg_permutation' / run_id / 'asm/contigs.fasta.to_GENE.cs.paf'
    assert cs.is_file() and Path(str(cs) + '.evidence.json').is_file(), 'Missing completed cs PAF: ' + str(cs)

print('\nRunID\tTargetPairsDesign\tTargetPairsFastq\tActual/Design')
for r in rows(wd / 'result/rrmg_afum_target_reads_genome_coverage_summary.tsv'):
    try:
        designed, actual = float(r['TargetPairsDesign']), float(r['TargetPairsFastq'])
        ratio = actual / designed if designed > 0 else float('nan')
        print(f"{r['RunID']}\t{designed:g}\t{actual:g}\t{ratio:.4f}")
        assert actual == designed, f"Target pair mismatch for {r['RunID']}: {actual} != {designed}"
    except (KeyError, ValueError, TypeError):
        print('[WARN] Missing/non-numeric target read counts:', r['RunID'])
        raise

for plan_row in plan:
    run_id = plan_row['run_id']
    run_dir = wd / 'temp/rrmg_permutation' / run_id
    with (run_dir / 'design/read_mix_audit.tsv').open(newline='') as fh:
        audit = {r['key']: r['value'] for r in csv.DictReader(fh, delimiter='\t')}
    total = int(round(float(plan_row['total_pairs'])))
    background = int(audit['background_pairs_written'])
    target = int(audit['target_pairs_written'])
    assert int(audit['total_pairs_written']) == total
    assert background + target == total
    assert int(audit['gzip_members_per_output']) == 1
    assert len(audit['mix_generation_id']) == 32
    fastp_json = run_dir / 'qc' / f"{plan_row['sample']}.Seed{plan_row['seed']}.fastp.json"
    report = json.loads(fastp_json.read_text())
    observed_reads = int(report['summary']['before_filtering']['total_reads'])
    assert observed_reads == 2 * total, f'fastp input mismatch for {run_id}: {observed_reads} != {2 * total}'
print('[OK] Planned sample IDs and required output files are present.')
print('[OK] Every target FASTQ count, mixed-read audit and fastp input count matches the plan.')
PY

# 主要汇总文件：
#   $WD/result/rrmg_target_gene_map_merged.tsv
#   $WD/result/rrmg_reads_to_gene_depth_summary.tsv
#   $WD/result/rrmg_target_only_reads_to_gene_depth_summary.tsv
#   $WD/result/rrmg_afum_target_reads_genome_coverage_summary.tsv
# 深度曲线还需要每个 run/readmap/*.depth.tsv；coverage 窗口图还需要
# target_readmap/*.bins.tsv。它们已保留在测试目录中，没有删除。
# VarN=0 是现有流程占位，不代表已完成基因突变检测。
# 本流程测试不等于真实检出性能评估；审计通过只证明数据流完整。
# 可以继续生成下一段的 QC 预览，帮助查看结果，但不能据此替代正式重复设计。


# ===== 13. 汇总全部绘图输入到 WD/figure_data（服务器或 Ubuntu） =====
# 本段不运行 R、不重新比对、不读取 FASTQ；只复用 Step3/4/5 已完成的 TSV、PAF 和 depth。
# figure_data 是可搬运的绘图输入快照；figure_output 留给本地 R 生成，不混放分析输入。
# 四张 result 汇总表会先重新收集，再按完整 plan 严格核对。任一样品缺失即停止，
# 不会把缺失值伪装成零。旧的非本工具管理 figure_data 会被保护，需先人工移走。
SITE_ARGS=()
if [[ -n "${SITES_SOURCE:-}" || -n "${CATALOGUE_SOURCE:-}" ]]; then
    test -s "$SITES_SOURCE"
    test -s "$CATALOGUE_SOURCE"
    SITE_ARGS=(--sites "$SITES_SOURCE" --catalogue "$CATALOGUE_SOURCE")
fi

bash "$CODE/rrmg_collect_figure_data.sh" \
    --wd "$WD" \
    --plan input/rrmg/rrmg_permutation_plan.tsv \
    --root temp/rrmg_permutation \
    --gene "$WD/ref/CYP51A_wt.fasta" \
    --require-base-evidence \
    --out "$WD/figure_data" \
    "${SITE_ARGS[@]}"

test -s "$WD/figure_data/COMPLETE.json"
test -s "$WD/figure_data/FILE_MANIFEST.tsv"
python3 -m json.tool "$WD/figure_data/COMPLETE.json"
find "$WD/figure_data" -type f | wc -l
du -sh "$WD/figure_data"

# figure_data 包含：完整 plan、三张汇总表、每个样品的 contig PAF、两类 reads 的
# PAF/depth、target-genome bins、目标基因 FASTA、可选位点注释、校验清单和 COMPLETE.json。
# contig 无命中时 PAF 可以是合法空文件；其状态由汇总表和 manifest 区分，不按缺失处理。
# 可选打包后下载：
#   cd "$WD" && tar -czf rrmg_figure_data.tar.gz figure_data


# ===== 14. 下载 figure_data 后，在本地 Windows PowerShell 运行全部 R 绘图 =====
# 这一段是 PowerShell 示例，不要粘贴到 Linux Bash。
# 先通过 WinSCP/scp/rsync 将服务器的整个 figure_data 下载到本地工作目录，例如：
#   D:\RRMG\SpikeIN\figure_data
# 同时在 Windows clone/download 本仓库。随后在 PowerShell 运行：
#
# $CODE = 'D:\RRMG\rrmg-fungal-spikein'
# $WD = 'D:\RRMG\SpikeIN'
# & "$CODE\rrmg_plot_all_windows.ps1" -WorkDir $WD -CodeDir $CODE
#
# 多个 R 版本并存时显式指定：
# & "$CODE\rrmg_plot_all_windows.ps1" -WorkDir $WD -CodeDir $CODE `
#     -Rscript 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe'
#
# 输出固定写入：
#   $WD\figure_output\local_preview     总体检出、breadth/depth和设计比较
#   $WD\figure_output\alignment_sites   contig/read blocks和位点覆盖
# 两个 R 程序读取同一个 figure_data 快照，但生成不同图组；不会运行任何生信比对。


# ===== 15. 失败后如何检查/续跑（参考命令，不自动执行） =====
# 日志：$WD/logs/local/<本次调用>/summary.tsv 和 row_1.log / row_2.log。
# 软件详细日志在对应 run 的 reads/、qc/、asm/、readmap/ 目录。
# 下列示例仅在确认需要时取消注释；不要从头重新生成已有 plan。
#
# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate fungamr_meta
# CODE="/absolute/path/to/rrmg-fungal-spikein"
# export RRMG_CONFIG="$HOME/rrmg_local_test_01/rrmg_config.sh"
# source "$RRMG_CONFIG"
# bash "$CODE/rrmg_submit.sh" step3 2 -t 4 -m 16 --jobs 1
#
# 默认只复用带当前 read_mix/fastp/target-pair 审计且 generation ID 一致的输出。
# 旧结果缺少 read_mix_audit.tsv 时，先用新版 Step1 给现有背景计数，再对该 row
# 依次强制重跑 Step2、Step3、Step4/5。重做 reads 后不能混用旧 assembly/mapping。
