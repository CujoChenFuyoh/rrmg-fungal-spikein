# 安装、从零准备输入与跨服务器运行

版本 0.3.0。样品目录、实验 plan 和分析主体保持原结构；本次改变的是调度入口和软件路径解析。

**重要变化：`rrmg_submit.sh` 不再递交 PBS，而是直接运行任务。**
激活环境后运行 `./rrmg_submit.sh step2 1-100`，默认顺序处理原 plan 中的 run_index 1–100。
不要在 HPC 登录节点直接执行重分析；把同一句命令放进自己的 PBS 作业。
本地入口不再执行 `.pbs` 文件，也不再生成假的 PBS 环境变量。

## 1. 安装和配置

```bash
conda activate fungamr_meta
git clone https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git
cd rrmg-fungal-spikein
chmod +x rrmg_submit.sh rrmg_run_local.sh
python -m pip install .
rrmg --version
export CODE="$PWD"
export WD="$HOME/SpikeIN_local_test"
rrmg init --workdir "$WD" --pipeline-dir "$CODE"
export RRMG_CONFIG="$WD/rrmg_config.sh"
nano "$RRMG_CONFIG"
mkdir -p "$WD/input/rrmg" "$WD/ref" "$WD/temp" "$WD/result"
rrmg doctor --config "$RRMG_CONFIG"
```

已有克隆目录时不要覆盖重克隆。已有配置时不要再次 init，直接指定 `RRMG_CONFIG`。
pip 只安装流程代码，不安装分析工具，也不复制实验输入。也可从 Git URL 安装：

```bash
python -m pip install 'git+https://github.com/CujoChenFuyoh/rrmg-fungal-spikein.git@main'
rrmg init --workdir /path/to/SpikeIN
```

不指定 `--pipeline-dir` 时，核心代码冻结到 `.rrmg/core-0.3.0/`，防止以后的 pip 升级悄悄替换当前实验代码。
Git URL 只能安装已经推送的代码；本地改动使用 `pip install .`。
可用 `rrmg export /path/to/new_core` 导出带校验清单的核心代码和 examples；不会覆盖已有目录。
离线安装可先 `python -m pip wheel --no-deps . -w dist`，复制 wheel 后用 `pip install --no-index --no-deps 文件.whl`。
正式分析固定已推送的 commit SHA 和测试过的 Conda explicit lock。

## 2. 一个环境，路径由使用环境决定

`fungamr_meta` 内已有 datasets、Java 和分析软件时，不另建下载环境。下载器直接调用 PATH 中的 `datasets`，
Step3 也直接调用当前环境中的 fastp、SPAdes 和 minimap2，不再反复 `conda run`。
配置里的 Java/BBMap 路径默认留空；发行版不固定任何电脑的路径。

```bash
command -v datasets
datasets version
command -v java
readlink -f "$(command -v java)"
java -version
command -v randomreadsmg.sh
```

如果 BBMap 是另行安装、没有进入 PATH，可在私有配置中填写：

```bash
RANDOMREADSMG_SH="/absolute/path/bbmap/randomreadsmg.sh"
```

BBMap 目录需包含配套的 Java class 和辅助脚本；不能只复制一个 randomreadsmg.sh。
部分 PBS 系统需要在作业脚本中显式指定 Java。这个例外写在
[普通 PBS](../examples/pbs.pbs) 和 [PBS array](../examples/pbs_array.pbs) 的注释示例里，
放在 `conda activate` **之后**，按真实环境解析 `bin/java` 的符号链接；通常无需启用。
没有把任何机构特定路径写进发行版。

默认配置为小规模试跑：一般每任务 4 线程，Step2 为固定 seed 可重复性默认 1 线程，模拟 Java heap 4g、SPAdes 16 GB。
真实宏基因组可能需要更多资源，按作业申请调整。`--jobs 4` 会叠加四份资源，不是共用四线程。
doctor 只做工具/基础路径检查，不证明样品数据完整；尚未生成的 plan/reference 会提示 warning。

## 3. Prepare manifest → 下载 genome

下面只用于**新的独立试跑目录**。已有实验不要重新生成 manifest/plan 或更改 ID/随机种子。
示例背景表是两个物种的演示数据，不代表原始研究群落。实际分析改成你的物种丰度表。
输入列为 `Name`、`NCBI ID`、`NCBI Rank`、`Summed base count`；
可用 prepare 脚本的 `--help` 查看自定义列名参数。

```bash
source "$RRMG_CONFIG"
python3 "$PIPELINE_DIR/rrmg_prepare_manifest.py" \
  --input "$PIPELINE_DIR/examples/taxa_demo.tsv" \
  --outdir "$WD/input/rrmg" --topN 2 \
  --target-name 'Aspergillus fumigatus' --target-taxid 330879 \
  --target-genome-id Afum --target-fasta "$WD/ref/Afum_Af293.fasta"

bash "$PIPELINE_DIR/rrmg_download_genomes.sh" \
  --manifest "$WD/input/rrmg/community_manifest.tsv" \
  --out-manifest "$WD/input/rrmg/community_manifest.downloaded.tsv"

cat "$WD/input/rrmg/download_failed.tsv"
test -s "$WD/ref/Afum_Af293.fasta"
```

manifest 的六列是 `genome_ID name taxid rel_abund fasta role`，实际以 Tab 分隔。
目标行先为 `role=target, rel_abund=0`，spike 比例在生成 plan 时指定。
新示例用 [Af293 菌株 taxid 330879](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=330879)；
不能仅凭物种层级 taxid 就保证下载的是指定菌株。已有实验仍保留原 ID，单独核对参考序列身份。
下载器复用已有非空 FASTA；优先 reference，否则按组装级别回退。回退包可能包含很多 genome，
不是保证每次仅下载一个。`--threads` 目前不实现并行下载。
失败诊断在 `download_logs/`；背景下载失败会剔除并重新归一化，必须先审核失败列表。
目标下载失败会报错。正式分析保留 `--keep-tmp` 的 assembly report、FASTA 和 SHA-256。
目标基因 `ref/CYP51A_wt.fasta` 必须另行提供并验证，本下载器不自动提取 CYP51A。

## 4. 测试 genome 与生成小规模 plan

```bash
source "$PIPELINE_DIR/rrmg_environment.sh"
bash "$PIPELINE_DIR/rrmg_test_randomreadsmg_genomes.sh" \
  --manifest "$WD/input/rrmg/community_manifest.downloaded.tsv" \
  --out-manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
  -w "$WD" --test-pairs 200 --readlen 150 --seed 999

python3 "$PIPELINE_DIR/rrmg_make_plan.py" \
  --manifest "$WD/input/rrmg/community_manifest.tested.tsv" \
  --out "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  --sample-prefix Pilot_Afum \
  --target-name 'Aspergillus fumigatus' --target-taxid 330879 \
  --target-genome-id Afum --target-fasta "$WD/ref/Afum_Af293.fasta" \
  --target-abundances '0.01' --total-pairs '1000000' --seeds '1-2' --readlen 150

python3 "$PIPELINE_DIR/rrmg_make_background_plan.py" \
  --plan "$WD/input/rrmg/rrmg_permutation_plan.tsv" \
  -w "$WD" --outdir "$WD/temp/rrmg_background" --background-seed 1000
```

结果是 1 个背景条件、2 个 permutation。兼容性测试会剔除失败背景并归一化；
先核对 `randomreadsmg_genome_test/` 下的 summary/log，再决定是否接受群落变化。
低深度用于检查能否执行，不要求一定检出目标基因。扩大规模前核对 FASTQ 实际对数、工具版本和参考基因组。
本更新没有修正历史 reads 数量异常，也没有改变随机种子和 abundance 语义。
部分 RandomReadsMG 版本把 `reads=N` 解释为每个 contig 的 reads 数，而不是整个 genome 的总数。
多 contig 目标尤其要核对 Step5 的 `TargetPairsFastq` 与 `TargetPairsDesign`，不能直接把设计丰度当作实际丰度。

## 5. 相同入口用于直接运行、普通 PBS 和 array

```bash
conda activate fungamr_meta
cd "$PIPELINE_DIR"
./rrmg_submit.sh step1
./rrmg_submit.sh step2 1-2
./rrmg_submit.sh step3 1-2
./rrmg_submit.sh step4 1-2
./rrmg_submit.sh step5 1-2
# 或在 Step3 后用 step45 1-2 一起做 Step4/5。
# 可选：./rrmg_submit.sh step3b 1-2
```

step1 的数字是 background plan 的 `bg_index`，其他步骤是原 permutation plan 的 `run_index`。
不是 Seed、目标物种 ID，也不是包含表头的文件行号。单数字只跑一个 ID，省略则跑表中全部 ID。
支持保留原编号的子集表；不自动重新编号。任何缺失/重复 ID 在开始计算前报错。
`--dry-run` 只验证选择范围和入口，不会检查所有 FASTQ 是否完整。
`rrmg_run_local.sh` 和 `rrmg run` 仍可用；三者共用同一核心分发逻辑。

普通 PBS 示例的计算部分是：

```bash
conda activate fungamr_meta
./rrmg_submit.sh step2 1-100 --jobs 1
```

支持 array 时，每个子任务只运行自己的 ID：

```bash
ROW="${PBS_ARRAY_INDEX:-${PBS_ARRAYID:-}}"
./rrmg_submit.sh step2 "$ROW" --jobs 1
```

**不能写成 `1-$ROW`**，否则每个 array 子任务都会重复处理前面样品。
完整头部、环境初始化、可选 Java 设置和提交命令见 [examples/README.md](../examples/README.md)。
PBS 的 queue/account/resource 语法由用户按服务器调整，核心流程没有固定 PBS 要求。
旧根目录 `_array.pbs` 只保留为兼容示例；新本地入口不再依赖它们。

## 6. 继续旧实验、汇总和备份

```bash
conda activate fungamr_meta
export RRMG_CONFIG=/path/to/existing/rrmg_config.sh
rrmg run step2 409-480 --dry-run
rrmg run step2 409-480 --jobs 1
rrmg run step2 225-230 --jobs 1
# 相应完成后，按同样范围运行 step3，再 step45。
```

原输入、背景 FASTQ/index、plan 中的绝对路径仍需有效，clone/pip 不会迁移这些数据。
旧配置的 `/EDIT_ME/` 或失效 Java 路径需清空或修正；在调用前主动激活环境。
保留 `STEP2_KEEP_TARGET_READS=1` 供 Step5 使用。不要同时操作同一个样品，或重建仍被消费的背景。
默认沿用各步骤原有复用条件；`--force` 才显式重跑，可能删除指定样品已有 assembly/binning。
每次执行都有 `logs/local/.../summary.tsv` 和 `row_N.log`，部分失败时继续尝试其余选定样品，最终返回非零。
Step3 是 fastp QC，不是 FastQC；Step5 的 coverage bins 不是 VAMB genome bins。
Step3 的空 VCF/`VarN=0` 仍是占位，不能当作完成 fungAMR 突变检测。

绘图输入收集命令见 [README 的 Collect 部分](../README.md#7-resume-collect-and-preserve-results)。
无需全部重复结束才能收集，但必须对照 plan 核对缺失 ID，不能把缺失当零。
只做 reads/contigs 时不用等 VAMB；按图需要额外备份 per-base depth 和 genome coverage window 表。
不要在任务运行期间升级代码或更改配置。
