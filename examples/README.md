# 运行示例：调度与分析分开

核心 `rrmg_submit.sh` 执行指定步骤与原始任务ID；环境由调用方激活，
调度资源由用户自己的作业脚本设置。

## 直接运行

```bash
conda activate fungamr_meta
export RRMG_CONFIG=/path/to/SpikeIN/rrmg_config.sh
bash ./rrmg_submit.sh step2 1-100 -t 1 -m 4 --jobs 1
# 等价的示例包装：bash examples/local.sh step2 1-100
```

## PBS 与 PBS array

两种通用模板集中在 [pbs_example/](../pbs_example/README.md)：

- [普通 PBS](../pbs_example/pbs.pbs)：按队列申请资源，在作业内执行一个ID或范围。
- [PBS array](../pbs_example/pbs_array.pbs)：每个task执行一个原始ID。

两者都先读取 `$WD/rrmg_config.sh`，再初始化Conda和环境。
初始化命令、PBS头部、提交方法和可选Java设置见该目录的README。
