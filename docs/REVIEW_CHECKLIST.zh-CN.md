# 仓库精简与重复计算确认清单

核对基线：`9776fdcd9942769808d133c4a4f076a51e983b4f`，版本0.3.1。
当前版本0.3.4；保留原Git提交历史，不更改当前PolyU部署文件、计划或既有实验结果。

## 已实施

- [x] README收束为入口和导航；安装细节、输出定义、绘图和故障各有唯一详细说明。
- [x] 新增快速上手：配置一次、原任务ID、每步小规模资源、收集与绘图路径。
- [x] runbook移除个人E/F/H盘和用户R库默认路径；采用用户参数/环境检测。
- [x] pilot默认单样品并发、4线程/16GB组装；不默认 `--force`。
- [x] 背景规划程序直接输出scheduler-free命令，移除教程中的“忽略旧qsub提示”。
- [x] R注释描述当前图层与颜色含义；移除修改任务口吻和过期手工验收数字。
- [x] 保留两个通用PBS模板；按确认移除根目录五个旧stage array wrappers。
- [x] 移除旧target-only绘图补丁；统一使用当前alignment准备及COMPLETE输出契约。
- [x] 拒绝重复seed、重复设计、执行表中重复run_id及不一致的copy_number模式。
- [x] VAMB日志不再依赖未导出的Conda环境名。
- [x] 按本项目的完整绘图目标保留PAF；不新增统一 `--no-paf` 主入口，不关闭默认PAF。
- [x] 保留read-pair计数、fastp输入审计、generation、坐标校验、缺失值及输出目录保护。

## 重复入口：是否必要

| 项目 | 核对结果 | 处理 |
| --- | --- | --- |
| `rrmg run`、`rrmg_submit.sh`、`rrmg_run_local.sh` | 共用实现；是安装入口、兼容别名和实现层，不是三套算法 | 新用户只讲 `rrmg run`，兼容别名保留 |
| 普通PBS与array模板 | 两种调度模型，必要 | 保留 `pbs_example/` 两模板 |
| 根目录五个旧array wrappers | stage固定的历史便利入口，功能重叠 | 已移除；使用两种通用模板 |
| Step4、Step5与Step45 | 单独选择和组合便利；选择单步时另一分支被跳过 | 保留，不是重复执行两遍 |
| 旧target-only补丁与alignment重算 | Step4现已直接产生两类read mapping | 已移除常规重算；统一collector只验证和复制现有输出 |
| 两个R程序 | 总体统计和逐样品alignment展示，输出不同 | 均保留 |

## 重复比对：是否必要

| 比对 | 核对结果 | 待确认动作 |
| --- | --- | --- |
| Step4 mixed-QC→gene BAM及PAF | 扫描同一对reads两次，分别供depth和blocks | 保留输出；以后验证单次比对同时导出或严格缓存 |
| Step4与绘图准备中的mixed-QC→gene | 旧版重复扫描 | 已消除；collector直接复用Step4 PAF/depth |
| Step3与绘图中的contig→gene | 旧版重新比对且口径不同 | 已复用Step3 PAF，并在绘图副本使用同一MAPQ>=20过滤 |
| target-only→gene | 是read-block图的独立来源 | 已纳入Step4并独立续跑 |
| Step5 target-only→genome | 参考不同，回答全基因组覆盖问题 | 可选，但不是上项的重复/替代 |

完整基因绘图常规路径为5次minimap2调用（Step3一次、Step4两类reads各两次；
不含Step5、可选species mapping/VAMB）。collector不调用比对程序。

## 仍需确认 / 下一轮实现

- [ ] 优先统一各阶段的参考、工具版本和参数指纹；修复改参考后旧完成标记仍被复用。
- [ ] 实现“保留assembly，只重做gene mapping”，不让改gene参考触发不必要的重新组装。
- [x] alignment绘图直接复用Step3 contig PAF，并采用与Step3汇总相同的MAPQ>=20过滤。
- [ ] 扩展参考、工具版本和参数指纹；当前PAF/depth复用仍以generation及结构校验为主。
- [ ] 拆出独立QC步骤，支持从零开始的纯read-based模式。
- [ ] 将 `VarN=0` 占位改成 `NA/not_performed`，同步collector和下游字段兼容。
- [ ] 明确研究的gene-hit判定阈值与对照设计；新增零丰度设计前验证空target-read路径。
- [x] 用户确认后移除五个旧PBS入口和旧target-only补丁。
- [ ] 作者选择许可证；当前不生成未经确认的LICENSE。
- [x] 用户指定本轮Step4/绘图收集整理版本为0.3.3；发布不表示其余待确认研究口径已改变。

## 验证范围

本轮测试记录应绑定基线提交和修改文件SHA-256；本地验证记录在交付清单中。
不沿用旧文档的“29项/39项测试”“50 included”等数字作为当前通过证明。
测试使用临时小型输入和mock工具；真实HPC、完整规模及未执行的R渲染须单独说明。
