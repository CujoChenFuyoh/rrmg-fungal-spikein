# 故障与续跑

## 命令没有提交PBS

这是当前入口的设计：`rrmg run` / `rrmg_submit.sh` 在当前环境执行。
需要排队时，用自己站点的外壳提交；参考 `pbs_example/`。
0.3.2已移除根目录五个旧 `_array.pbs`；请使用通用模板调用 `rrmg run`。

## 计划被拒绝

- 使用原 `run_index`；不要把seed、行号或新编号当作任务ID。
- `1-100` 要求其中每个ID均在plan中。存在编号缺口时，按现有ID分段/逐行选择。
- 重复seed、重复设计或重复 `run_id` 会共享输出目录，必须先修正输入，不能靠并发运行区分。
- `copy_number` 的计划与执行尚未统一，当前明确拒绝；使用经核对的 `read_fraction`。

## fastp输入计数不符或旧输出被拒绝

旧版压缩流直接拼接可能使fastp提前停止。当前Step2重写单一gzip member，
并记录精确pair数和generation；Step3/4检查fastp实际输入计数。
没有审计的旧结果应先保存。确认目标ID后，按Step1背景审计、Step2重建、
Step3和所需Step4/5重建的顺序处理，不能把新reads与旧QC/assembly拼接。
`--force` 可能删除该阶段已有输出，尤其Step3的assembly；先备份再显式使用。

## 改了参考或参数，流程仍跳过

当前完成标记未覆盖所有参考/参数指纹。这是已登记的限制，不代表新参考已经运行。
固定现有实验的参考/参数；需要比较新设置时使用独立项目或核对后显式重建。
“保留assembly、只重比对”的独立入口仍在确认清单中，当前不能靠修改config推断完成。

## Windows FASTA 导致contig覆盖百分比偏低

0.3.2的Step3长度计算已忽略FASTA空白字符；LF和CRLF输入采用相同的实际碱基数。
此修复不改变contig命中、PAF区间、MAPQ规则或组装结果。仅安装更新不会修正旧汇总表：
已有结果需要保留reads、QC和assembly，重新生成contig比对汇总、收集表及受影响的统计图。
不要直接用 `rrmg run step3 --force` 代替仅重比对，因为它会重建assembly。

## 位点文件报目录、空路径或坐标错误

`--sites`、`--catalogue` 应成对指定为真实非空TSV文件；不标注位点时同时省略。
用 `test -s "$FILE"` 核对变量和文件。图形准备会核对gene FASTA方向、长度、
codon坐标与参考碱基，并输出排除原因；不使用氨基酸位置乘3推算含内含子的坐标。

## 图没有命中或旧图仍在

先区分输入缺失和真实零覆盖：查看汇总状态、逐碱基depth与PAF，而非仅看图层颜色。
绘图输出采用固定目录和管理清单，只有本次成功完成才替换上一版。
已有 `COMPLETE.json` 可能属于上一轮；要以本次 `[DONE]` 为准。
目录中有用户修改/新增文件时先保存或移动目录，不删除保护清单绕过检查。

## Windows R

Linux原生R直接运行 `.R`；WSL调用Windows R时使用 `rrmg_windows_r.sh`。
可显式设置 `WIN_RSCRIPT` 和可选 `WIN_R_LIB`。不指定库时采用Windows R自己的库配置。
具体检查命令和依赖见绘图说明。不要让Linux R加载Windows编译的包。
