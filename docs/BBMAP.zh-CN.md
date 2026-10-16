# BBMap/BBTools：安装、下载与路径核对

RRMG 需要 BBTools 的 `randomreadsmg.sh`（背景）和 `randomreads.sh`（目标）。
`git clone` / `pip install rrmg-fungal-spikein` 只复制流程代码，不安装这些程序。
`environment.yml` 已列出 `bbmap`，但已有的 `fungamr_meta` 不一定由这份文件创建。
先检查实际环境，再选择下面一种安装方式；不需要另建 BBMap 环境。

## 1. 判断正在调用哪套 BBMap

在 Linux/WSL/服务器中执行：

```bash
conda activate fungamr_meta
printf 'Environment: %s\n' "$CONDA_PREFIX"
conda list -n fungamr_meta bbmap
command -v randomreadsmg.sh || true
command -v randomreads.sh || true
command -v java
```

`conda list` 的包名列表里有 `bbmap` 才说明它由该环境的 Conda 管理。
`command -v` 找到程序不等于它安装在环境里：也可能是全局 PATH 或独立目录。
Bioconda 的 launcher 可以是符号链接，因此还应检查 `readlink -f` 的真实位置。

还需核对实际项目配置，因为显式 override 优先于 PATH：

```bash
export RRMG_CONFIG=/absolute/path/to/SpikeIN/rrmg_config.sh
source "$RRMG_CONFIG"
source "$PIPELINE_DIR/rrmg_environment.sh"
printf 'Background tool: %s\nTarget tool: %s\nJava: %s\n' \
    "$RANDOMREADSMG_SH" "$RANDOMREADS_SH" "$JAVA_BIN"
readlink -f "$RANDOMREADSMG_SH"
readlink -f "$RANDOMREADS_SH"
```

`RANDOMREADSMG_SH` / `RANDOMREADS_SH` 若已经指定，优先使用它们；否则
`BBMAP_DIR` 会先加入 PATH，再查找工具。Java 的解析是独立的，可以来自
`fungamr_meta`，同时让 BBMap 本身位于环境外。

## 2. 方式 A：安装进现有 fungamr_meta

[Bioconda 提供 BBMap 包](https://bioconda.github.io/recipes/bbmap/README.html)。
以下版本是本机独立安装所报告的版本，用于迁移示例，不是“最新版本”要求；
若你原服务器使用其他版本，请先核对并替换，不要在续跑时无记录地升级。

```bash
conda activate fungamr_meta
BBMAP_VERSION=39.33

# 先审阅求解结果，不改变已安装软件。
conda install -n fungamr_meta --override-channels \
    -c conda-forge -c bioconda --strict-channel-priority \
    "bbmap=$BBMAP_VERSION" --dry-run

# 确认不会不希望地更新/降级已有依赖后，再执行安装。
conda install -n fungamr_meta --override-channels \
    -c conda-forge -c bioconda --strict-channel-priority \
    "bbmap=$BBMAP_VERSION"

conda list -n fungamr_meta bbmap
command -v randomreadsmg.sh
command -v randomreads.sh
readlink -f "$(command -v randomreadsmg.sh)"
readlink -f "$(command -v randomreads.sh)"
java -version
```

本仓库没有验证所有服务器上该版本与现有依赖的 Conda 求解结果。
发生冲突时先检查 dry-run，不要为了安装而直接重建正在使用的环境。
`environment.yml` 是新环境的起点，不是对已有环境运行自动更新的指令。

如果项目原先绑定了另一套 BBMap，在实际使用的 `rrmg_config.sh` 中把以下
三项清空，重新开一个终端并激活环境，然后重新核对真实路径：

```bash
BBMAP_DIR=""
RANDOMREADSMG_SH=""
RANDOMREADS_SH=""
```

不要以为装入 Conda 后就自动替换了项目的显式路径。确认解析到环境的
`bin/` 或 `share/` 后，再运行 `rrmg doctor --config "$RRMG_CONFIG"`。

## 3. 方式 B：官方归档独立下载，继续使用环境内 Java

当不希望 Conda 求解改变现有依赖时，可以下载完整 BBTools 归档。
[官方仓库](https://github.com/bbushnell/BBTools)指向 SourceForge 下载；
这里使用[固定 39.33 归档](https://sourceforge.net/projects/bbmap/files/BBMap_39.33.tar.gz/download)，
不使用会随时间改变的 `latest/download`。

在允许联网的登录节点/本地 Linux 下载；这不是模拟或 assembly 任务：

```bash
conda activate fungamr_meta
BBMAP_VERSION=39.33
SOFTWARE_ROOT="$HOME/software"
RELEASE_ROOT="$SOFTWARE_ROOT/BBMap-$BBMAP_VERSION"
ARCHIVE="$SOFTWARE_ROOT/BBMap_$BBMAP_VERSION.tar.gz"
mkdir -p "$SOFTWARE_ROOT"

# 防止覆盖以前的安装；目录已存在时检查/复用它，不直接重复解压。
test ! -e "$RELEASE_ROOT"
test ! -e "$ARCHIVE"
curl -fL --retry 3 \
    "https://sourceforge.net/projects/bbmap/files/BBMap_$BBMAP_VERSION.tar.gz/download" \
    -o "$ARCHIVE"

# 避免把失败下载的 HTML 页面当作软件包。
tar -tzf "$ARCHIVE" >/dev/null
sha256sum "$ARCHIVE" > "$ARCHIVE.sha256"
mkdir "$RELEASE_ROOT"
tar -xzf "$ARCHIVE" -C "$RELEASE_ROOT"
export BBMAP_DIR="$RELEASE_ROOT/bbmap"

test -s "$BBMAP_DIR/randomreadsmg.sh"
test -s "$BBMAP_DIR/randomreads.sh"
test -s "$BBMAP_DIR/current/shared/Shared.class"
bash "$BBMAP_DIR/bbversion.sh"
java -version
```

如果下载失败留下半个归档，先检查该文件，移走失败副本或使用新的下载目录后重试。
归档解压结构应有 `bbmap/`；以上文件检查失败时不要继续运行流程。
SHA-256 用于保存下载/传输一致性记录，不表示已经和独立公布的签名或摘要比对。

然后把路径写入私有 `rrmg_config.sh`；写在激活 `fungamr_meta` 后生效：

```bash
BBMAP_DIR="$HOME/software/BBMap-39.33/bbmap"
RANDOMREADSMG_SH="$BBMAP_DIR/randomreadsmg.sh"
RANDOMREADS_SH="$BBMAP_DIR/randomreads.sh"
# JAVA_BIN / JAVA_HOME 通常留空，从激活环境的 PATH 解析 Java。
```

不能只复制两个 `.sh` 文件。需要完整目录，包括 Java classes、配套 shell
脚本及资源。另一种等价的离线迁移方式是复制已经验证过的整个 BBMap 目录，
在新服务器更新这三个路径，并重新检查版本/Java/文件摘要。

不建议为了让 `conda list` 出现 BBMap 而额外安装第二份；环境外安装本身有效。
对续跑来说，选择正确且明确的一套程序比安装方式更重要。

## 4. 迁移与 PBS 检查

独立目录不会被 Conda explicit lock 包含，因此额外保存版本和 launcher/class
摘要；仅有相同版本号不证明整个目录内容相同：

```bash
source "$RRMG_CONFIG"
source "$PIPELINE_DIR/rrmg_environment.sh"
BBMAP_ACTUAL_DIR="$(dirname "$(readlink -f "$RANDOMREADSMG_SH")")"
bash "$BBMAP_ACTUAL_DIR/bbversion.sh"
sha256sum "$RANDOMREADSMG_SH" "$RANDOMREADS_SH" \
    "$BBMAP_ACTUAL_DIR/current/shared/Shared.class"
conda list -n fungamr_meta --explicit
rrmg doctor --config "$RRMG_CONFIG"
```

在计算节点的 PBS 作业内也做同样的工具路径检查，不要只看登录节点。
在解析工具和运行步骤前必须先 `conda activate fungamr_meta`；使用
`CONDA_PREFIX` 的配置项也应在激活后求值。没有另外的 `ncbi_datasets` 或
BBMap 环境。部分 PBS 作业需要显式 Java 路径，参考
[PBS 示例](../examples/README.md)，不固定为任何机构的特殊路径。

先完成少量任意 DNA 的兼容性测试，再运行正式样品。`doctor` 能检查基础路径，
但仅有 launcher/Java 存在不证明模拟结果或配对计数正确。
