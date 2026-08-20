# dualporo_HM XML 输入文件分类检查

检查时间：2026-07-16
远程目录：`/data/datafile/lsl/project/dualporo_HM`

## 结论摘要

该目录中共检查到 71 个 XML 输入文件。结合 XML 内容、运行输出、后处理脚本、分析数据清单和论文草稿判断，论文主结果真正使用的输入文件是 `5.15hyst` 系列中的 F/R 两组工况：

- F 系列：改变 WAG 周期/slug size，水/CO2 质量比约为 1。
- R 系列：固定 WAG 周期为 3 个月，改变水/CO2 质量比。

其余 XML 大多是逐步建立模型、调试求解器、测试注入量/井位/单孔隙对照/短期运行的探索文件。它们不建议作为论文主结果证据使用，但可以保留作模型演化、敏感性测试或排错记录。

## 判断依据

主要依据如下：

1. `artical/drafts/论文初稿.md` 明确写到：本文仅将 F 系列和 R 系列纳入结果分析，其他探索注入体积或井位变化的算例不作为本文证据。
2. `_postproc_5.15series/plot_final.py` 中最终绘图列表只包含：
   - `FREQ = 5.15hyst_F1_1mo, 5.15hyst, 5.15hyst_F6_6mo, 5.15hyst_F12_12mo, 5.15hyst_Fcont_5yr`
   - `RATIO = 5.15hyst_R0p25, 5.15hyst_R0p5, 5.15hyst, 5.15hyst_R2, 5.15hyst_R3`
3. `storage_compare_article_F_series` 和 `storage_compare_article_R_series` 是论文图表和数据来源，里面的 `analysis_manifest.json`、`near10_summary_metrics.csv`、`injection_summary.csv` 与上述 F/R 工况一致。
4. 最终主算例 XML 均包含双重介质、裂缝连续介质、WAG `SourceFlux`、滞回相渗表和对应 PVD/run 输出。
5. 早期 `0.x`、`1.x`、`4.12`、`4.13`、`5.14` 等目录显示出明显的逐步建模和调参痕迹，且与论文最终 F/R 数据清单不一致。

## 最终使用的主案例 XML

### F 系列：WAG 周期/slug size

这些是论文中 WAG 周期对比的主输入文件。除基准 `5.15hyst` 外，文件名已经直接标出周期。

| 论文工况 | 含义 | XML 文件 | 关键设置 | 判断 |
|---|---|---|---|---|
| F1 | 1 个月 WAG 周期 | `5.15hyst_F1_1mo/input.xml` | CO2 `-20.0 mol/s`，水 `-48.89 mol/s`，水/CO2 质量比约 1，`maxTime=630720000` | 最终使用 |
| F3 / R1 | 3 个月基准 | `5.15hyst/input.xml` | CO2 `-20.0 mol/s`，水 `-48.89 mol/s`，水/CO2 质量比约 1，`maxTime=315360000` | 最终使用，F/R 共用基准 |
| F6 | 6 个月 WAG 周期 | `5.15hyst_F6_6mo/input.xml` | CO2 `-20.0 mol/s`，水 `-48.89 mol/s`，水/CO2 质量比约 1，`maxTime=630720000` | 最终使用 |
| F12 | 12 个月 WAG 周期 | `5.15hyst_F12_12mo/input.xml` | CO2 `-20.0 mol/s`，水 `-48.89 mol/s`，水/CO2 质量比约 1，`maxTime=630720000` | 最终使用 |
| F60 / Fcont | 60 个月，即 5 年半周期 | `5.15hyst_Fcont_5yr/input.xml` | CO2 `-20.0 mol/s`，水 `-48.89 mol/s`，水/CO2 质量比约 1，`maxTime=630720000` | 最终使用 |

说明：`storage_compare_frequency_F1_F3_F6_F12_Fcont_near10/README.md` 也确认了 F1、F3、F6、F12、Fcont 这五个案例，其中 F3 对应 `5.15hyst`。

### R 系列：水/CO2 质量比

这些是论文中水/CO2 质量比对比的主输入文件。WAG 周期固定为 3 个月。

| 论文工况 | 水/CO2 质量比 | XML 文件 | 水注入通量 | 判断 |
|---|---:|---|---:|---|
| R0.25 | 0.25 | `5.15hyst_R0p25/input.xml` | `-12.2225 mol/s` | 最终使用 |
| R0.5 | 0.5 | `5.15hyst_R0p5/input.xml` | `-24.445 mol/s` | 最终使用 |
| R1 / F3 | 1.0 | `5.15hyst/input.xml` | `-48.89 mol/s` | 最终使用，F/R 共用基准 |
| R2 | 2.0 | `5.15hyst_R2/input.xml` | `-97.78 mol/s` | 最终使用 |
| R3 | 3.0 | `5.15hyst_R3/input.xml` | `-146.67 mol/s` | 最终使用 |

`storage_compare_article_R_series/injection_summary.csv` 中的 CO2 通量固定为 `20.0 mol/s`，水通量按上述比例变化，且所有工况 CO2 注入开启总时间为 `157680000 s`。

## 与最终结果配套的后处理和论文文件

这些不是输入 XML，但用于确认哪些案例真正进入论文：

| 路径 | 用途 |
|---|---|
| `artical/drafts/论文初稿.md` | 论文正文草稿，明确限定主结果为 F/R 两组工况 |
| `artical/drafts/论文_6页版_模板重排.pdf` | 排版版论文，可作为论文参考版本 |
| `artical/data/analysis_manifest.json` | F/R 结果指标清单 |
| `storage_compare_article_F_series/` | F 系列论文图表和 CSV 数据 |
| `storage_compare_article_R_series/` | R 系列论文图表和 CSV 数据 |
| `_postproc_5.15series/extract_series.py` | 从 VTK/PVD 提取双重介质 CO2 状态、压力、断层应力指标 |
| `_postproc_5.15series/plot_final.py` | 最终图表脚本，硬编码 F/R 最终案例列表 |

注意：`论文初稿.md` 对断层滑移指标表述较保守，强调当前论文主要采用孔隙压力指标；`论文_6页版_模板重排.pdf` 摘要中则写入了断层面应力投影结果。输入文件分类不受这个差异影响，但后续写论文时需要统一安全性表述。

## 过时或探索性 XML 分组

下面这些文件或目录不建议作为论文主结果使用。

### 0.x 初始模型搭建

| 目录 | 判断 | 原因 |
|---|---|---|
| `0.1basemodel/` | 过时 | 早期 base model，通量极小，未形成最终 WAG/F/R 设计 |
| `0.2basemodel16core/` | 过时 | 主要用于核心数量/并行测试 |
| `0.3basemodelmaxDt/` | 过时 | 主要用于最大时间步调试 |
| `0.4gasinjectionAquiferShortTimeBase/` | 过时 | 早期气体注入测试 |
| `0.5gasAlterInjectionAInquifer/` | 过时 | 早期气水交替可行性探索，未进入论文数据链 |

这些目录通常缺少最终 `run.log` 证据，XML 中也没有最终的滞回相渗和 F/R 参数体系。

### 1.x 断层、重力、partial dual 建模阶段

| 目录 | 判断 | 原因 |
|---|---|---|
| `1.1fault_alter_400*200_base/` | 过时 | 断层 400x200 基础模型建立 |
| `1.2fault_alter_400*200_base_gravity/` | 过时 | 加入重力和分区 |
| `1.3fault_alter_400*200_base_gravity_partialDual/` | 过时 | partial dual 阶段 |

这些是最终几何和双重介质模型的前置演化步骤，但不是论文最终对比工况。

### 4.12 / 4.13 / 5.14 中间调参阶段

| 目录 | 判断 | 原因 |
|---|---|---|
| `4.12wag_for_storage_status/` | 被 5.15 取代 | 早期 WAG 统计基准，CO2/水通量为 `-200/-488.89 mol/s`，比最终 5.15 主案例高 10 倍 |
| `4.12F1_freq1mo/`, `4.12F6_freq6mo/`, `4.12F12_freq12mo`, `4.12Fcont_noWAG/` | 被 5.15 取代 | 早期频率系列，无最终滞回设置，未进入论文最终数据 |
| `4.12R0p5_wg0p5/`, `4.12R2_wg2/`, `4.12R3_wg3/` | 被 5.15 取代 | 早期气水比系列，通量尺度与最终不同 |
| `4.12Vscale_*` | 调参 | 注入量/孔隙度/求解器时间步调试 |
| `4.13/`, `4.130.5t/` | 中间版本 | 4.12 后继续调参，但不在论文 F/R 数据清单中 |
| `5.14lowKfhighPhim/` | 中间版本 | 低裂缝渗透率/高基质孔隙度尝试，随后被 5.15hyst 系列取代 |

这些目录可用于追溯模型如何收敛到最终方案，但不应混入论文主结果。

### 5.15hyst 派生敏感性和调试案例

| 目录 | 判断 | 原因 |
|---|---|---|
| `5.15hyst_0p01PV_*` | 探索性 | 0.01 PV 注入量/周期测试，不在论文 F/R 主结果 |
| `5.15hyst_0p05PV_*` | 探索性 | 0.05 PV 注入量/周期测试，不在论文 F/R 主结果 |
| `5.15hyst_0p05PV_topprod_*` | 探索性 | 顶部生产井/生产位置变体 |
| `5.15hyst_amg/` | 求解器测试 | AMG 配置测试，PVD 很短，不能作为主结果 |
| `5.15hyst_iluk1/`, `5.15hyst_iluk2/` | 求解器测试 | ILU(k) 线性求解器配置测试 |
| `5.15hyst_fv003_auto/`, `5.15hyst_fv003_forced/` | 参数调试 | fracture volume 或相关参数测试 |
| `5.15hyst_fperm_e12/`, `5.15hyst_fperm_e13/` | 敏感性 | 裂缝渗透率量级测试，晚于论文主数据生成 |

这些可以作为附加敏感性分析或排错材料，但不属于论文主线案例。

### 单孔隙、短时、低通量和井位变体

| 目录 | 判断 | 原因 |
|---|---|---|
| `5.15single_*` | 对照/探索 | 单孔隙模型，不是论文主张的双重介质主结果 |
| `t2_*` | 短时或替代井位测试 | `maxTime=63072000 s`，约 2 年，不是论文近 10 年 F/R 主结果 |
| `injbot_*` | 注入井位置/周期/短时测试 | 井位或短时变体，其中 `injbot_R*_2yr` 为 2 年 |
| `q0p1_*`, `q3_*`, `dpq0p1_*`, `testtime/` | 低/高通量和时间测试 | 主要用于通量、时间长度或排错，不在论文数据清单 |
| `2d_*`, `2dq01_*` | 输出结果目录/派生测试 | 没有直接作为最终 F/R XML 主输入列入论文 |

## 建议保留方式

建议将目录逻辑上分成三类：

1. `final_for_paper`：保留 `5.15hyst`、`5.15hyst_F*`、`5.15hyst_R*` 和 `storage_compare_article_*`。
2. `supporting_analysis`：保留 `artical/`、`_postproc_5.15series/`、`storage_compare_frequency_*`。
3. `archive_or_exploration`：其余 `0.x`、`1.x`、`4.x`、调参和短时测试目录。

如果后续要清理目录，建议先不要删除旧算例，而是移动到 `archive/` 并保留这个报告，避免后处理脚本或论文图表中的相对路径失效。
