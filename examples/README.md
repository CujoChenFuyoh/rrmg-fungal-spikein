# 运行示例：调度与分析分开

这些模板需要先填写路径、资源、queue/account，不能原样提交。
核心 `rrmg_submit.sh` 只执行，不会 qsub，也不会激活/切换环境。

## 直接运行

```bash
conda activate fungamr_meta
export RRMG_CONFIG=/path/to/SpikeIN/rrmg_config.sh
bash ./rrmg_submit.sh step2 1-100 --jobs 1
# 等价的示例包装：bash examples/local.sh step2 1-100
```

## 普通 PBS（不需要 array）

修改 [pbs.pbs](pbs.pbs)，保留服务器验证过的 PBS 头部。
模板的 PBS Pro `select=1:ncpus=4:mem=24gb` 可按站点改为 Torque 所需的
`nodes=1:ppn=4`、`mem=24gb` 等。线程/内存必须与配置匹配，给 Java/SPAdes 留额外内存。
100 个样品顺序处理可能超过 walltime；首次用一个样品估计时间后再分块。

```bash
qsub -v RRMG_STAGE=step2,RRMG_ROWS=1-100 examples/pbs.pbs
# 下一阶段只在上游完成后提交：
# qsub -v RRMG_STAGE=step3,RRMG_ROWS=1-100 examples/pbs.pbs
```

## PBS array（仅支持 array 的服务器）

修改 [pbs_array.pbs](pbs_array.pbs)。每个 task 用一个原始 ID，不能用 `1-$ROW`。
PBS Pro 示例：

```bash
qsub -J 1-100 -v RRMG_STAGE=step2 examples/pbs_array.pbs
```

Torque 可能使用 `qsub -t 1-100` 和 `PBS_ARRAYID`，见
[官方 qsub 说明](https://docs.adaptivecomputing.com/torque/6-1-1/adminGuide/Content/topics/torque/commands/qsub.htm)。
模板兼容两种变量，但 array 语法、限流、依赖和队列限制按具体站点核对。
Step1 用背景表 bg_index，其他步骤用 permutation plan run_index；勿把 array 数量当成每个 task 的循环范围。

## 可选：在作业中显式选择 Java

多数情况只需 `conda activate fungamr_meta`。部分 PBS 系统的非交互作业环境需要显式指定 Java。
下面是环境内 Java 的参考写法，已作为注释放进两个 PBS 文件，不是发行版默认路径：

```bash
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate fungamr_meta
export JAVA_BIN="$(readlink -f "$CONDA_PREFIX/bin/java")"
export JAVA_HOME="$(dirname "$(dirname "$JAVA_BIN")")"
export PATH="$JAVA_HOME/bin:$PATH"
java -version
```

Conda 的 `bin/java` 可能是指向 `lib/jvm/bin/java` 的符号链接，因此先解析真实路径。
也可以明确指定另一套兼容的 Java，但需在计算节点验证 `java -version` 和小规模模拟。
使用上述 export 时，私有配置里的 Java 选项应留空/继承环境，不能再次写入过期路径。
datasets 始终从同一个已激活环境调用，不需要单独的 ncbi_datasets 环境。

