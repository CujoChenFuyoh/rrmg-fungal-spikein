#!/usr/bin/env bash
# RRMG 0.3.1 本地 Ubuntu/WSL 逐段运行手册
# ------------------------------------------------------------
# 用法：在编辑器打开本文件，按编号依次复制到同一个 Ubuntu Bash 终端。
# 带反斜杠的多行命令、if/fi、Python heredoc 必须作为一个完整块执行。
# 不在 Windows PowerShell 里执行，不需要 qsub，不调用 PBS 文件。
# 建议首次逐段检查，不直接 bash runbook.sh 一次跑到底。
# 本文件是本机手册，开头路径不是发行版默认设置。
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
CODE="/mnt/e/github_repo/rrmg-fungal-spikein"
export WD="$HOME/rrmg_local_test_01"
BBMAP_DIR_LOCAL="/mnt/f/DB/bbmap"
# 上面是本机独立 BBMap 目录，不是 fungamr_meta 中的软件安装路径。
# 复制到服务器时不能照搬；按 docs/BBMAP.zh-CN.md 安装/下载并填写真实目录。
# 若使用 Conda 内 BBMap，先激活环境再按下面的命令查找目录；不要在激活前解析 PATH。
# BBMAP_DIR_LOCAL="$(dirname "$(readlink -f "$(command -v randomreadsmg.sh)")")"
GENE_SOURCE="/mnt/e/genedata/spikeIn/ref/CYP51A_wt.fasta"
# runbook、R 绘图脚本和辅助脚本均随仓库提供，统一通过上面的 CODE 调用。
SITES_SOURCE="/mnt/f/DB/FungAMR_catalogue/CYP51A_FungAMR_sites.tsv"
CATALOGUE_SOURCE="/mnt/f/DB/FungAMR_catalogue/Afum_fungamr.catalogue.tsv"
export WIN_RSCRIPT="/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
export WIN_R_LIB="C:/Users/q6321/AppData/Local/R/win-library/4.6"
export RRMG_CONFIG="$WD/rrmg_config.sh"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate fungamr_meta

# 报错即停。不要把报错后的命令也粘贴进去。
set -eo pipefail
cd "$CODE"
test "$(tr -d '\r\n' < VERSION)" = "0.3.1"
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
# 使用本地代码，因为本地 0.3.1 修改可能尚未推送到 GitHub。
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
rrmg export "$WD/code_snapshot_0.3.1"


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
# 如背景规划脚本末尾仍打印旧 qsub 提示，忽略它；按下面的新入口运行。
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
bash "$CODE/rrmg_submit.sh" step1 1 -t 4 -m 4 --jobs 1 --force
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
bash "$CODE/rrmg_submit.sh" step2 1-2 -t 1 -m 4 --jobs 1 --force
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
#   STEP3_SKIP_QC=0             # 设 1 仅在已有可验证 QC reads 时跳过 fastp
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
bash "$CODE/rrmg_submit.sh" step3 1-2 -t 12 -m 100 --jobs 2 --force
find "$WD/temp/rrmg_permutation" -path '*/asm/contigs.fasta' -type f -size +0c -print
# 应看到两个 contigs.fasta。无目标基因命中不等于 assembly 运行失败。
# 重新拼接后继续第 11–15 段，重收集 contig 汇总并重建绘图输入；不能只重画旧汇总。


# ===== 9. Step4：QC reads → CYP51A，得到 per-base depth/coverage =====
# 设置参数（线程直接写在下面命令；无需编辑 config）：
#   -t 4 / --threads 4          # 每个样品的 minimap2/samtools 线程数
#   GENE_REF="ref/CYP51A_wt.fasta"  # 相对 WD 或绝对路径；必须与后续绘图参考一致
# Step4 没有可靠的程序总内存参数，因此不接受 -m；PBS mem 在作业脚本中按服务器设置。
# 提交范围和 --jobs 规则与 Step2 相同；PBS ncpus 至少为 jobs × -t。
# Step4 依赖 Step3 的 QC reads；--force 只重建本步骤 mapping/depth 输出。
bash "$CODE/rrmg_submit.sh" step4 1-2 -t 12 --jobs 2 --force


# ===== 10. Step5：target-only reads → Af293 genome coverage =====
# 设置参数（线程直接写在下面命令；无需编辑 config）：
#   -t 4 / --threads 4                # 每个样品的 minimap2/samtools 线程数
#   TARGET_GENOME="ref/Afum_Af293.fasta"  # 相对 WD 或绝对路径
#   STEP45_MAKE_BINS=1                   # 生成固定窗口 coverage 表
#   STEP45_BIN_SIZE=10000                # coverage 窗口大小（bp），不是 MAG bin
# Step5 没有可靠的程序总内存参数，因此不接受 -m；PBS mem 在作业脚本中按服务器设置。
# Step5 必须有 STEP2_KEEP_TARGET_READS=1 生成的 target-only 双端 reads。
# 提交范围和 --jobs 规则与 Step2 相同；提高 jobs 时 CPU、内存和磁盘 I/O 同时增加。
bash "$CODE/rrmg_submit.sh" step5 1-2 -t 12 --jobs 2 --force
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
             'rrmg_afum_target_reads_genome_coverage_summary.tsv']:
    table = rows(wd / 'result' / name)
    ids = [r['RunID'] for r in table]
    print(f'{name}: {len(ids)} rows; missing={sorted(expected-set(ids))}; extra={sorted(set(ids)-expected)}')
    assert set(ids) == expected and len(ids) == len(expected), 'Incomplete/duplicate summary: ' + name
assert not missing_files, 'Missing output files:\n' + '\n'.join(missing_files)

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
#   $WD/result/rrmg_afum_target_reads_genome_coverage_summary.tsv
# 深度曲线还需要每个 run/readmap/*.depth.tsv；coverage 窗口图还需要
# target_readmap/*.bins.tsv。它们已保留在测试目录中，没有删除。
# VarN=0 是现有流程占位，不代表已完成基因突变检测。
# 本流程测试不等于真实检出性能评估；审计通过只证明数据流完整。
# 可以继续生成下一段的 QC 预览，帮助查看结果，但不能据此替代正式重复设计。


# ===== 13. 生成可直接查看的 figure output：PNG + PDF + 绘图源数据 =====
# 分析已完成时，只需恢复第 0 段变量，运行第 11–15 段；不必重跑 Step1–5。
# 在 Ubuntu/WSL 中调用 Windows 的 Rscript.exe，使用 Windows 的 R 包库。
# 不调用 /usr/bin/Rscript，也不让 Linux R 加载 Windows 编译的包。
# 两边路径由包装脚本自动转换；检查依赖，不自动安装/修改 R 包。
# 指标/分组参考原 Windows 版 RRMG_permutation_summary_plots.R 和
# RRMG_Afum_target_reads_genome_coverage_plots.R；图形用 base R 重绘，不依赖其旧路径。
bash "$CODE/rrmg_windows_r.sh" --check
test -s "$CODE/rrmg_plot_preview.R"

# 固定目录。新图全部生成成功后替换工具管理的旧输出，过期分页不会混入新结果。
# 若目录里有手动添加/修改的文件，脚本会拒绝覆盖；请先保存或移动该目录。
FIG_OUT="$WD/figure_output/local_preview"
bash "$CODE/rrmg_windows_r.sh" "$CODE/rrmg_plot_preview.R" "$WD" "$FIG_OUT"

printf '\nFigure output: %s\n' "$FIG_OUT"
ls -lh "$FIG_OUT"/*.png "$FIG_OUT"/*.pdf
column -t -s $'\t' "$FIG_OUT/sample_status.tsv"
cat "$FIG_OUT/warnings.txt"

# 主要图片（样品/设计条件多时自动分页）：
#   overview_01.png / .pdf：gene breadth、gene mean depth、target genome breadth、target pairs 审核。
#   gene_depth_01.png / .pdf：各重复的目标基因逐碱基深度曲线。
#   design_hit_rates.png / .pdf：contig/read gene-hit 比例热图，显示有效数/计划数。
#   contig_metrics_by_design_01.*：contig breadth、covered bp、alignment length、identity。
#   read_metrics_by_design_01.*：mapped alignments、breadth、mean/max depth、breadth >=5x。
#   target_genome_metrics_by_design_01.*：target genome 预期/实测深度、breadth、深度比。
#   target_expected_vs_observed.*：target genome 预期/实测深度及 target pairs 预期/实际散点图。
#   summary_source_data.tsv / gene_depth_source_data.tsv：图片对应的实际数值。
#   detection_by_design.tsv：每个设计的 hits、有效分母、计划数与缺失数。
#   sample_status.tsv / warnings.txt：缺失、target pairs 不一致、深度表与汇总不一致。
# NA 不补成 0；真正测得的 0 会保留。不同 reference 序列不拼成一条曲线。
# 检出比例仅是当前样品的描述值，分母是可评价结果数；不是总体灵敏度。
# 每个条件的 n 从实际表格计算；不足 3 个有效值时只画散点，不画箱线。
# 不绘制 VarN 占位值，不把未运行的 Kraken/VAMB 当作阴性结果。
# 下方第 14–15 段补齐原 alignment-block/FungAMR-site 图的输入并出图。
# 注意 target-only reads -> gene 不等于 Step5 的 target-only reads -> genome。
# 若需要更多样品，改变 plan 后重新收集并运行本段；默认绘制 plan 中的所有 RunID。
# 部分样品缺失时仍可作 QC 预览，但必须查看状态表，不能视为完整实验。
#
# Ubuntu 桌面打开概览图（可选，取消注释后运行）：
# xdg-open "$FIG_OUT/overview_01.png"
# WSL 在 Windows 文件资源管理器打开结果目录（可选）：
# explorer.exe "$(wslpath -w "$FIG_OUT")"
#
# 同一个 R 脚本也能在 Windows 使用，下面是 PowerShell 示例，不要在 Bash 粘贴：
# $Rscript = 'C:/Program Files/R/R-4.6.0/bin/Rscript.exe'
# $CODE = 'E:/github_repo/rrmg-fungal-spikein'
# $WD = 'H:/genedata/SpikeIN'
# $FIG_OUT = "$WD/figure_output/local_preview"
# & $Rscript "$CODE/rrmg_plot_preview.R" $WD $FIG_OUT
# Invoke-Item "$FIG_OUT/design_hit_rates.png"
# 支持 result/ 三张表及 figure_data/{contig_mapping,read_to_gene,target_genome}/ 备份结构。
# 默认以 input/rrmg/rrmg_permutation_plan.tsv 为完整待检查队列；若该文件不存在，
# 才用 figure_data/metadata/selected_plan.tsv，且在 warnings.txt 说明仅为选定子集。


# ===== 14. 补齐位点标注和 alignment-block 绘图输入 =====
# 仍在 Ubuntu 的 fungamr_meta 环境运行，不提交 PBS，不重新合成 reads/assembly。
# 必须保留 Step2 target-only reads 和 Step3 QC reads；此步骤不会补造已删除的数据。
# 使用既有 contigs、QC reads、target-only reads 重新对同一 gene FASTA 做比对，
# 保证三类 PAF/depth 对应同一参考序列。不覆盖已有 Step1–5 文件。
# 本段重复执行会重新生成所选 rows 的整套绘图输入，不是追加；仅重新出图可跳到第 15 段。
# 位点表来自已存在的 TSV，并用 catalogue 的 chromosome/codon_positions 和
# 当前 gene FASTA 的真实参考碱基核验；不使用 aa_pos*3，不推断或修复错误坐标。
# catalogue status 非 ok、参考密码子不符或坐标不符的条目进入排除审计表。
# 注释仅为参考位点的覆盖展示，不代表相应突变已检出；promoter/TR 不在 gene 区间时不补画。
# 恢复变量，避免单独复制本段时空字符串被 Python 当作当前目录 '.'。
SITES_SOURCE="${SITES_SOURCE:-/mnt/f/DB/FungAMR_catalogue/CYP51A_FungAMR_sites.tsv}"
CATALOGUE_SOURCE="${CATALOGUE_SOURCE:-/mnt/f/DB/FungAMR_catalogue/Afum_fungamr.catalogue.tsv}"
: "${CODE:?先设置 CODE 为仓库目录}"
: "${WD:?先设置 WD 为已有分析工作目录}"
for INPUT_FILE in "$WD/ref/CYP51A_wt.fasta" "$SITES_SOURCE" "$CATALOGUE_SOURCE"; do
    if [[ ! -f "$INPUT_FILE" || ! -s "$INPUT_FILE" ]]; then
        printf '[ERROR] 输入必须是非空文件：%s\n' "$INPUT_FILE" >&2
        exit 1
    fi
done
ALIGN_INPUT="$WD/figure_data/CYP51A_alignment"
python3 "$CODE/rrmg_prepare_alignment_inputs.py" \
    --wd "$WD" \
    --gene "$WD/ref/CYP51A_wt.fasta" \
    --out "$ALIGN_INPUT" \
    --rows 1-2 \
    --threads 4 \
    --sites "$SITES_SOURCE" \
    --catalogue "$CATALOGUE_SOURCE"

# --rows 是 plan 中的 run_index；扩大队列时可用 1-100 或 409-480,225-230，或 all。
# COMPLETE.json 只在所有所选样品完成后生成。无命中允许空 PAF，但 depth 必须逐碱基完整。
# 只有看到本次 [DONE] 后才继续。失败不会替换上一轮成功结果，不要将旧 COMPLETE 当作本次成功。
test -s "$ALIGN_INPUT/COMPLETE.json"
find "$ALIGN_INPUT" -maxdepth 2 -type f \( -name '*.paf' -o -name '*.paf.gz' -o -name '*.depth.tsv' \) -print
column -t -s $'\t' "$ALIGN_INPUT/site_annotation_audit.tsv"
# 每个 row_XXXXX/：contig PAF、mixed-QC PAF.gz + depth、target-only PAF.gz + depth。
# 总目录：gene FASTA、plot_manifest.with_target_only.tsv、校验后的位点表、原始表与审计记录。
# BAM 只是此步骤的中间文件，验证 depth 后删除；原 BAM、FASTQ 和 contigs 均不删除。


# ===== 15. Windows R 生成位点标注、contig alignment、read-block 图 =====
# 使用旧脚本的 ggplot2/patchwork 图形结构与橙色/灰色 reads 配色，已移除旧 setwd。
# 用当前的 Windows R 包库，不需要在 Ubuntu 安装这些 R 包。
ALIGN_INPUT="$WD/figure_data/CYP51A_alignment"
ALIGN_FIG_OUT="$WD/figure_output/alignment_sites"
bash "$CODE/rrmg_windows_r.sh" \
    "$CODE/rrmg_plot_alignment.R" "$ALIGN_INPUT" "$ALIGN_FIG_OUT"

find "$ALIGN_FIG_OUT" -maxdepth 2 -type f \( -name '*.png' -o -name '*.pdf' -o -name '*.tsv' \) -print
# alignment_depth/：每个样品参考位点 + contig blocks + 两类 reads depth，PNG/PDF 和多页 PDF。
# read_blocks/：每样品两类 reads 的 alignment blocks，PNG/PDF、多页 PDF 和 PAF 统计表。
# CYP51A_FungAMR_site_coverage_by_run.tsv / ..._by_design.tsv：每位点覆盖统计，非突变检出表。
# read blocks 默认最多显示每来源 350 条，是显示抽样；位点深度仍用完整 depth 数据。
# 可选打开：explorer.exe "$(wslpath -w "$ALIGN_FIG_OUT")"


# ===== 16. 失败后如何检查/续跑（参考命令，不自动执行） =====
# 日志：$WD/logs/local/<本次调用>/summary.tsv 和 row_1.log / row_2.log。
# 软件详细日志在对应 run 的 reads/、qc/、asm/、readmap/ 目录。
# 下列示例仅在确认需要时取消注释；不要从头重新生成已有 plan。
#
# source "$HOME/miniconda3/etc/profile.d/conda.sh"
# conda activate fungamr_meta
# CODE="/mnt/e/github_repo/rrmg-fungal-spikein"
# export RRMG_CONFIG="$HOME/rrmg_local_test_01/rrmg_config.sh"
# source "$RRMG_CONFIG"
# bash "$CODE/rrmg_submit.sh" step3 2 -t 4 -m 16 --jobs 1
#
# 默认只复用带当前 read_mix/fastp/target-pair 审计且 generation ID 一致的输出。
# 旧结果缺少 read_mix_audit.tsv 时，先用新版 Step1 给现有背景计数，再对该 row
# 依次强制重跑 Step2、Step3、Step4/5。重做 reads 后不能混用旧 assembly/mapping。
