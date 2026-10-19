# 快速上手

本页提供最小运行路径。输入下载和参数细节见[安装与输入教程](INSTALL.zh-CN.md)，
逐段试跑见[runbook.sh](../runbook.sh)。分析入口统一使用 `rrmg run`。

## 1. 安装与生成配置

在 Linux/Ubuntu/WSL 的 Bash 中执行；已有分析环境时直接激活，不重复安装软件。

```bash
conda activate fungamr_meta
git clone https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git
cd rrmg-fungal-spikein
python -m pip install .
export CODE="$PWD"
export WD="$HOME/rrmg_pilot"
rrmg init --workdir "$WD" --pipeline-dir "$CODE"
export RRMG_CONFIG="$WD/rrmg_config.sh"
rrmg doctor
```

已有 clone 时进入已有目录；已有 config 时跳过 `init`，直接设置 `RRMG_CONFIG`。
`pip` 安装代码；Conda 环境、基因组、FASTQ 和参考基因不会由它复制。

## 2. 核对输入

新实验按安装教程第3–4节依次完成 manifest、下载、兼容性测试、plan 和背景设计。
已经准备好输入的实验可直接核对以下文件：

| 输入 | 来源 |
| --- | --- |
| `input/rrmg/rrmg_permutation_plan.tsv` | `rrmg_make_plan.py` |
| 背景 `background_plan.tsv` / `background_index.tsv` | `rrmg_make_background_plan.py`；实际目录由config指定 |
| 背景和目标基因组 FASTA | manifest / plan 指向的真实文件 |
| `ref/CYP51A_wt.fasta` | 用户提供并核验的目标基因参考；config可指定其他名称 |

`run_index` 是原任务ID，不是 Seed 或含表头的文本行号。用原ID续跑，不能重新编号。
目前支持 `read_fraction`；`copy_number` 会明确报错。重复seed或重复输出目录也会报错。

## 3. 先跑两个重复

以下命令假定plan包含ID1、2。先检查选择范围，再执行：

```bash
rrmg run step2 1-2 --dry-run
rrmg run step1 -t 4 -m 4 --jobs 1
rrmg run step2 1-2 -t 1 -m 4 --jobs 1
rrmg run step3 1-2 -t 4 -m 16 --jobs 1
rrmg run step4 1-2 -t 4 --jobs 1
```

Step2合成和混合reads；Step3进行fastp、组装和contig mapping；Step4产生reads
覆盖度、depth及PAF。默认保留PAF，供read/contig alignment-block图使用。
Step3默认 `--only-assembler`，保留fastp，跳过BayesHammer。

每个阶段有 `logs/local/<调用目录>/summary.tsv` 和 `row_<ID>.log`。
只在上一阶段成功后继续；`--dry-run`检查选择和入口，不代替FASTQ完整性检查。
示例不默认使用 `--force`；需要覆盖时只对已核对的ID显式添加。

### 参数只改这几项

| 参数 | 含义 |
| --- | --- |
| `1-2` | 原任务ID范围；省略时运行表中全部ID，范围内缺失ID会报错 |
| `-t N` | 每个样品的线程数；Step2建议1以保持固定seed模拟行为 |
| `-m GB` | Step1/2 Java heap或Step3 SPAdes程序内存；无流程人为上限 |
| `--jobs N` | 当前进程内并发样品数；CPU/内存需求按并发数增加 |

## 4. 收集和绘图

```bash
source "$RRMG_CONFIG"
python3 "$PIPELINE_DIR/rrmg_collect_map_results.py" \
  --summary "$WD/result/rrmg_split_pipeline_summary.tsv" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --root "$WD/temp/rrmg_permutation" \
  --out "$WD/result/rrmg_target_gene_map_merged.tsv" --wd "$WD"
python3 "$PIPELINE_DIR/rrmg_collect_reads_to_gene_depth.py" \
  --root "$WD/temp/rrmg_permutation" \
  --out "$WD/result/rrmg_reads_to_gene_depth_summary.tsv"
```

随后按[绘图说明](../RUNBOOK_FIGURES_README.txt)生成总体图及alignment-block图。
后者需要contig→gene、mixed-QC→gene、target-only→gene三类PAF与对应depth；
图形输入准备完成后，仅重画图不需要再次运行比对。

Step5是另一个可选分析：`rrmg run step5 1-2 -t 4 --jobs 1`，
测量target-only reads对**整个基因组**的覆盖，不能替代target-only→gene的PAF。

## 5. 放到服务器

在站点PBS外壳中申请资源，读取同一config、激活Conda，再运行上面的 `rrmg run`。
模板见[pbs_example](../pbs_example/README.md)。一个普通PBS作业可以处理一个ID范围；
`--jobs`不提交PBS子任务。不要在HPC登录节点直接运行重计算。

当前PolyU续跑的固定Step3名单属于部署层，不由通用入口自动推断。
移植已有实验时保留plan、数据路径和这些选择表，避免重复组装已完成样品。
