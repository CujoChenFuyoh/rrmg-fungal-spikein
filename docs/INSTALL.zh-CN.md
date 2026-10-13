# 安装与迁移：继续使用原 RRMG 流程

0.2.0 只增加安装和命令行入口。`rrmg run` 实际执行原来的
`rrmg_run_local.sh`；样品目录、分析步骤和 plan 不变。不需要之前的续跑补充包。

## 一、git clone 后安装

下面的 `<YOUR-REPOSITORY-URL>` 必须换成你真实的 GitHub 仓库地址。本地改动尚未上传时，远程旧版本不会自动获得这些功能。

```bash
conda activate fungamr_meta
git clone <YOUR-REPOSITORY-URL> rrmg-fungal-spikein-github
cd rrmg-fungal-spikein-github
python -m pip install .
rrmg --version
```

核心代码安装到当前 Python/Conda 环境，新增 `rrmg` 命令。
**pip 不安装 BBMap、Java、SPAdes、minimap2、samtools 等外部分析程序**，也不复制参考序列、背景 reads、plan 或结果。已有分析环境直接沿用。

如果没有环境，先从克隆目录运行 `conda env create -f environment.yml`，再激活并执行 `python -m pip install .`。Conda YAML 不是精确版本锁；验证环境后保存 `conda list --explicit`。

## 二、直接从 GitHub 安装

```bash
# 把地址和提交 SHA 换成真实值；目标提交必须已包含 pyproject.toml。
python -m pip install "git+https://github.com/OWNER/REPOSITORY.git@COMMIT_SHA"
rrmg --version
```

可使用已发布的 tag，但完整 commit SHA 更明确。服务器需要 Git 和仓库访问权限。
本项目目前没有自动发布到 PyPI；不能据此假定 `pip install rrmg-fungal-spikein` 已可从公共索引下载。

无外网服务器可在可联网机器准备 wheel，再复制到服务器：

```bash
python -m pip wheel --no-deps . -w dist
# 复制 dist/rrmg_fungal_spikein-0.2.0-py3-none-any.whl 到服务器，然后：
python -m pip install --no-index --no-deps rrmg_fungal_spikein-0.2.0-py3-none-any.whl
```

这只是纯 Python 安装/调度代码的 wheel，不包含分析软件二进制。

## 三、已有工作目录和配置：不需要 init

例如 UBDA 已经有原数据和配置：

```bash
conda activate fungamr_meta
export RRMG_CONFIG=/ubda/home/23124109r/rrmg-fungal-spikein-github/rrmg_config.sh

# 只检查任务范围/入口；不运行分析。
rrmg run step2 451-455 --jobs 1 --dry-run

# 以下分析命令应在计算节点或普通 PBS 作业中运行。
rrmg run step2 451-455 --jobs 1
rrmg run step3 451-455 --jobs 1
rrmg run step45 451-455 --jobs 1
```

行号只是示例，换成实际未完成的 `run_index`。原入口要求保留完整、连续编号的 plan；不要剪出几行后重新编号。不是 `Seed`，也不是包含表头的文本行号。

也可以显式指定配置：

```bash
rrmg run step3 451 --config /path/to/rrmg_config.sh --jobs 1
```

配置优先级：`--config` → `RRMG_CONFIG` → `--workdir/rrmg_config.sh` → 当前目录的 `rrmg_config.sh`。
若设置了 `RRMG_CONFIG`，它会优先于 `--workdir`；切换项目时推荐显式 `--config`。
`run` 和 `doctor` 使用该配置的 `PIPELINE_DIR`，所以不会绕过你选定的原代码版本。

## 四、新环境：只初始化部署配置，不重建实验

先激活目标分析环境，让工具路径自动检测来自正确环境。

### 继续绑定现有 GitHub 克隆目录

```bash
rrmg init \
  --workdir /ubda/home/23124109r/SpikeIN \
  --pipeline-dir /ubda/home/23124109r/rrmg-fungal-spikein-github \
  --env fungamr_meta
```

只新增工作目录下的 `rrmg_config.sh`。原 plan、背景、reads 和 assembly 不变。
已存在配置时拒绝覆盖，直接使用已有配置；也可用 `--config /path/new_config.sh` 明确创建另一份配置。

### 没有单独的代码克隆目录（从 Git URL 或 wheel 安装）

```bash
rrmg init --workdir /path/to/SpikeIN --env fungamr_meta
```

会把**原始核心脚本**复制到工作目录的 `.rrmg/core-0.2.0/`，写配置指向它，并生成 `rrmg_core_manifest.json`。
这不是另一套分析流程，只是冻结当前脚本副本，防止以后 pip 升级悄悄替换当前项目的核心代码。
初始化中断留下的目录不会自动删除或覆盖；检查完整性后可用 `--pipeline-dir` 显式绑定。

自动检测不到的软件路径会写成 `/EDIT_ME/...` 并提示。检查配置中的全部路径和资源参数后再运行：

```bash
rrmg doctor --workdir /path/to/SpikeIN
```

可以在初始化时显式提供路径：

```bash
rrmg init --workdir /path/to/SpikeIN \
  --conda-sh /path/to/conda/etc/profile.d/conda.sh \
  --bbmap-dir /path/to/bbmap \
  --randomreadsmg /path/to/bbmap/randomreadsmg.sh \
  --java-bin /path/to/java/bin/java
```

`doctor` 仍调用原 preflight：它检查工具和基本文件，不代表 reads 完整、实验设计正确或所有分析已完成。
配置是可执行 Bash 文件，只加载你信任的配置。

## 五、需要一份可直接复制的核心脚本

```bash
rrmg export /path/to/new_rrmg_core
```

输出原来的 `.sh`、`.pbs`、`.py`、配置模板、环境文件、说明和 SHA-256 清单。
已有目标目录会被拒绝，避免混入不同版本的脚本；不提供批量覆盖选项。
这份导出是运行脚本集合，不是新的 Git 仓库或可再次 `pip install .` 的源码发行包；需要安装源码时使用 git clone 或源码压缩包。

导出以后仍可按旧方式运行：

```bash
export RRMG_CONFIG=/path/to/SpikeIN/rrmg_config.sh
bash /path/to/new_rrmg_core/rrmg_run_local.sh step2 451-455
```

相应配置的 `PIPELINE_DIR` 必须指向这一份代码。只移动代码目录，不会自动改写配置里的旧绝对路径。

## 六、不支持 array 的服务器

`rrmg run` 永远调用本地入口，不执行 qsub。可直接放在已有普通 PBS 脚本中：

```bash
# 放在已有普通 PBS 作业的计算部分；保留站点验证过的资源头。
source /path/to/conda/etc/profile.d/conda.sh
conda activate fungamr_meta
export RRMG_CONFIG=/path/to/rrmg_config.sh

for ROW in 451 452 453; do
    rrmg run step2 "$ROW" --jobs 1 || exit 1
    rrmg run step3 "$ROW" --jobs 1 || exit 1
    rrmg run step45 "$ROW" --jobs 1 || exit 1
done
```

原 `rrmg_submit.sh` 仍是 **array 提交器**，不适用于不支持 array 的服务器。
脚本名称包含 `_array.pbs` 不影响本地模式：本地入口用 Bash 执行它，PBS 注释头不会被调度器解析。
每个样品的线程/内存由配置控制，普通 PBS 请求必须匹配；增大 `--jobs` 会叠加资源需求。

## 七、规划与汇总工具

保留所有原命令，也可从安装包调用：

```bash
rrmg path
rrmg script rrmg_make_plan.py -- --help
rrmg script rrmg_collect_map_results.py -- --help

# 读取现有结果并输出原来的 result/ 汇总表：
rrmg script rrmg_collect_supplementary_downstream.sh -- /path/to/SpikeIN
```

`script` 使用当前安装包的工具，不自动加载项目配置或激活环境；工具参数和输出覆盖行为与原脚本相同。
对于已固定旧代码版本的项目，应按需直接调用该项目 `PIPELINE_DIR` 下的工具。
补充汇总脚本已修复代码目录与数据目录分开时找不到 Python collector 的问题。

## 八、迁移与升级边界

- pip/clone 解决代码分发，不能代替迁移 FASTA、背景 FASTQ、plan 和已有结果。
- 表格里的 `/rds/...` 等旧绝对路径仍需核对，不自动重建 plan 或修改随机种子。
- 原 Step3/Step45 的输入检查、结果复用规则和 `--force` 行为没有改变。
- 不假定全部已完成；收集器仍只汇总找到的结果，样品完整性需对照 plan。
- Step3 的 `VarN=0` 仍是占位，并非完成 fungAMR 突变检测；已有 reads 数量异常也没有在此安装更新中修正。
- 推荐固定代码提交和已验证的 Conda lock。不要在作业运行期间 `pip install --upgrade`、`git pull` 或修改配置。
- 默认快照不会自动升级；绑定 Git checkout 的项目会受到后续 checkout 修改影响。
- native Windows 可安装、导出、初始化、查看帮助；完整分析要在 Linux/WSL 上运行。

打包遵循 [setuptools 包内数据支持](https://setuptools.pypa.io/en/stable/userguide/datafiles.html) 和
[pip Git 安装方式](https://pip.pypa.io/en/stable/topics/vcs-support/)。项目没有新增许可证授权，公开分享前请由维护者自行选择合适许可证。
