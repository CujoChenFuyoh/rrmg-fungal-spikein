# PBS 参考：统一读取 rrmg_config.sh

本目录集中普通 PBS 和 array 两个模板。先按服务器要求修改 PBS 头部及模板内
`RRMG_CONFIG` 的绝对路径。PBS 会读取该文件，初始化 Conda、激活环境、进入
`WD`，再通过 `PIPELINE_DIR` 调用核心流程；不用 `-V` 继承提交终端环境。
PBS 可能将脚本复制到 spool 目录，因此不要用脚本自身位置推断 config 路径。

## 首次生成 config

安装 RRMG 后，用现有命令生成 `$WD/rrmg_config.sh`，并检查生成内容：

```bash
rrmg init --workdir /path/to/SpikeIN \
  --pipeline-dir /path/to/rrmg-fungal-spikein \
  --conda-sh /path/to/miniconda3/etc/profile.d/conda.sh --env fungamr_meta
# 单独安装 BBMap 时可增加 --bbmap-dir /path/to/bbmap
```

已有 config 不会被覆盖。`CONDA_SH`、`RRMG_CONDA_ENV`、`WD` 和 `PIPELINE_DIR`
只在这份 config 中维护；队列、核数和 walltime 在 PBS 头部设置。
通用 `rrmg init` 使用通用背景目录；PolyU 续跑包使用自己的初始化脚本，仍写入
相同文件名 `rrmg_config.sh`，以保留其背景目录和续跑选项。

## 普通 PBS：按队列申请资源，不需要 array

[pbs.pbs](pbs.pbs) 使用 `q2s01`、`nodes=1:ppn=24` 和168小时作为写法示例，
未固定 `mem=`。队列名、资源语法及实际内存以你的服务器分配为准。

```bash
qsub -v RRMG_STAGE=step2,RRMG_ROWS=1-100,RRMG_THREADS=1,RRMG_MEMORY_GB=4 pbs_example/pbs.pbs
```

一个作业顺序处理指定的原始 ID 范围。需要更多并行作业时，分别提交互不重叠的
范围，或用普通 `qsub` 循环逐样品提交。各阶段应在对应上游成功后执行。
Step3 可改为适合 assembly 的队列和核数，再设置匹配的 `RRMG_THREADS` 与
`RRMG_MEMORY_GB`。`-m` 只控制程序内存，不会替代调度器内存申请。

## 支持 array 的 PBS

[pbs_array.pbs](pbs_array.pbs) 给出 PBS Pro `select` 资源写法，按站点修改队列
及资源。PBS Pro 与 Torque 的提交语法分别为：

```bash
qsub -J 1-100 -v RRMG_STAGE=step2,RRMG_THREADS=1,RRMG_MEMORY_GB=4 pbs_example/pbs_array.pbs
# Torque，仅在服务器支持array时使用，并将PBS头部改成该站点资源语法：
# qsub -t 1-100 -v RRMG_STAGE=step2,RRMG_THREADS=1,RRMG_MEMORY_GB=4 pbs_example/pbs_array.pbs
```

一个 task 只处理 `PBS_ARRAY_INDEX` 或 `PBS_ARRAYID` 指定的一个原始 ID，
不会执行 `1-$ROW`。Step1 对应 `bg_index`，其他步骤对应 `run_index`。
队列限流、array范围和依赖语法需按服务器支持情况设置。

## 可选 Java 设置

一般使用激活环境中的 Java。部分 PBS 系统需要显式选择 Java，两个模板都保留了
激活后的注释示例；先解析 `CONDA_PREFIX/bin/java` 的真实路径，再设置
`JAVA_HOME`。也可以在唯一的 config 中指定 `JAVA_BIN`，供环境辅助脚本读取。
datasets 及其他程序仍使用同一个 `fungamr_meta` 环境。
