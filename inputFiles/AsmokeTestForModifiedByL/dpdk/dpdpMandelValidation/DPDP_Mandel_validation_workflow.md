# DPDP Mandel FIM/Seq 验证流程

本文档用于在双重介质流固耦合程序改动后，重新验证本目录中的 FIM 与 Sequential
Mandel 算例。目标是让后续更新验证算例时有固定流程，避免混用历史试验文件、覆盖旧结果，
或把阶段性调通结果误认为最终物理验证通过。

2026-07-16 清理后，目录中只保留两个当前主验证 XML 输入：

1. `DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml`：同源 pressure 标定的全耦合法输入。
2. `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml`：同源 pressure 标定的迭代耦合法输入。

历史试错 deck 已删除，差异和删除原因记录在
`DPDP_Mandel_trial_history_and_cleanup.md`。技术结论和问题状态记录在
`DPDP_Mandel_findings_current_status.md`。Seq/FIM 耦合实现流程见
`DPDP_Mandel_coupling_flow_Seq_vs_FIM.md`。物性、均匀化、storage、输入模式和 XML 参数定义见
`DPDP_Mandel_theory_input_reference.md`。

## 1. 验证基准

### 1.1 问题定义来源

`problem_description/` 保存原始问题定义资料：

- `Mehrabian和Abousleiman - 2014 - Generalized Biot's theory and Mandel's problem of multiple‐porosity and multiple‐permeability poroel-800571.pdf`
- `mandel_problem2014.md`

该目录只放论文原文、问题描述、边界条件和参数来源说明。不要把运行结果、试错记录或最新结论写入
`problem_description/`；这些内容分别写入 workflow、trial history 或 findings。

### 1.2 解析解基准

`analitical_result/fig5c_primary_analitical.csv` 是手动取得的 Fig. 5c primary/matrix
压力解析解点。

`analitical_result/fig5c_secondary_analitical.csv` 是手动取得的 Fig. 5c secondary/fracture
压力解析解点。

这两个 CSV 是本验证流程的压力主基准，不应当被脚本重新生成或覆盖。脚本解析解
`script/dpdp_mandel_analytical.py` 可作为辅助检查，但最终误差表应以这两个手动 CSV 为准。
当前 findings 记录了手动 Fig. 5c CSV 与解析脚本曲线不一致的问题；在该问题解决前，不能把
解析脚本生成的位移加载函数视为独立物理基准。

`analitical_result/fig5d_analitical.csv` 可作为应力辅助诊断。当前位移驱动算例中，应力需要由
有效应力和 Biot 压力项重构，不作为主要通过/失败判据。

`analitical_result/` 中只长期保留手动解析 CSV 和带时间戳的当前主对比图。旧阶段图片应清理到
`DPDP_Mandel_trial_history_and_cleanup.md` 的记录中，而不是继续散放在该目录。

### 1.3 网格

使用 10x10 网格作为标准验证网格。当前已经认为 10x10 满足该 Mandel 验证精度要求，无需使用
20x20 或 40x40 作为常规验证。

## 2. 推荐算例

### 2.1 FIM 同源 pressure 标定算例

基准文件：

`DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml`

验证前必须确认：

- `couplingType="FullyImplicit"`
- `enableFimCrossStorage="1"`
- `crossStorageOffDiagScale="1.0"`
- mesh 为 `nx="{ 10 }"`、`nz="{ 10 }"`
- 裂缝体积分数 `fractureVolumeFraction="0.03"`
- 裂缝渗透率采用体积加权后的 bulk flux contribution，即 `kappa_f = v_f * k_f`
- `defaultReferencePorosity` 输入 intrinsic 孔隙度；程序只将 `v_i` 乘入 storage/reference porosity。
  `permeabilityComponents` 不会再由 dual-continuum 初始化乘 `v_i`，因此本验证中必须直接输入 `kappa_i`
- 该 deck 使用同源 pressure 标定的初始压力和位移荷载，是当前主 FIM 验证输入
- `10000 s` 之后保持 `forceDt="300"`。`tau≈10^3-10^4` 的基质压力晚期快速下降段对时间步较敏感；
  粗步长会在对比图中形成非物理折线

### 2.2 Sequential 同源 pressure 标定算例

当前优先使用：

`DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml`

验证前必须确认：

- `couplingType="Sequential"`
- `useIntrinsicInput="0"`
- mesh 为 `nx="{ 10 }"`、`nz="{ 10 }"`
- `crossStorageOffDiagScale="1.0"`
- `subcycling="1"`
- `sequentialConvergenceCriterion="SolutionIncrements"`
- matrix 与 fracture 的 `ElasticIsotropic` K/G 均为 effective 输入值，当前为
  `K=4.514e8`、`G=3.108e8`
- matrix 与 fracture 的 `BiotPorosity` Biot 系数均为 effective 输入值，当前为
  `alpha_m=0.382`、`alpha_f=0.601`
- `DualContinuumCrossFlow` 中的 `intrinsicMatrixBiot`、`intrinsicMatrixBulkModulus`、
  `intrinsicFractureBiot`、`intrinsicFractureBulkModulus` 保留 intrinsic 物理参数，
  用于等效/储量公式检查，不应误改成 effective 参数
- `maxSequentialPressureChange="1.0e2"`。受控测试显示，修正 cross-storage offdiag split 后，
  `100 Pa` 外迭代容差可完整跑通且曲线指标基本不变；不要恢复到 `1.0e3`，该设置会导致部分
  时间步外耦合只迭代一次、下一步再补偿，在 `tau≈0.8-1` 附近形成压力锯齿
- 不要恢复旧的 `10 Pa` 外迭代容差和全程细时间步。该保守设置主要是旧 storage / cross-storage
  split 问题下的兜底；当前储量与 split 修正后，`100 Pa` 已能无切步跑完整程
- `10000 s` 之后保持 `forceDt="300"`。该设置用于消除 `tau≈10^3` 和 `tau≈10^4` 附近由 late-time
  时间步分段造成的基质压力折线；它不是 pressure relaxation，也不改变物理参数
- 该 deck 使用同源 pressure 标定的初始压力和位移荷载，是当前主 Seq 验证输入
- 孔隙度输入为 intrinsic 孔隙度，裂缝总体孔隙贡献按
  `fractureVolumeFraction * fracturePorosity` 理解和检查

## 3. 运行规范

### 3.1 输出目录

每次验证必须使用新的独立输出目录，避免覆盖历史结果。推荐格式：

```bash
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<case>_<date>_<tag>/
```

每个保留算例单独建子目录，例如：

```bash
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/fim_sameSourcePressure_mesh10_YYYYMMDD_01/
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/seq_sameSourcePressure_mesh10_YYYYMMDD_01/
```

每个子目录至少保留：

- 实际运行使用的 XML 副本
- GEOS log
- `pressure_matrix_history.hdf5`
- `pressure_fracture_history.hdf5`
- `displacement_history.hdf5`
- 若有应力诊断，则保留 `stress_history.hdf5`

### 3.2 运行命令模板

在仓库根目录运行。MPI 分解可按机器资源调整，但同一轮 FIM/Seq 对比应尽量保持一致。
当前 10x10 Seq 验证可以单核运行，便于排查和复现。

Seq 单核命令模板：

```bash
build-ubuntu-lsl-release/bin/geosx -x 1 -y 1 -z 1 \
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml \
  -o inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<seq_output_dir> \
  > inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<seq_output_dir>.log 2>&1
```

MPI 命令模板：

```bash
mpirun -np 4 build-ubuntu-lsl-release/bin/geosx -x 4 \
  -i <xml_file> \
  -o <output_dir> \
  > <output_dir>.log 2>&1
```

如需修改 XML 参数，应优先复制到 `/tmp` 或本次验证目录中形成临时 XML，不直接修改原始基准文件。

## 4. 对比流程

### 4.1 完整性检查

每个算例完成后先检查 log：

- 是否达到 `maxTime="308412.0"`
- 是否出现 nonlinear failure
- 是否出现过多 time step cut
- FIM Newton 迭代是否稳定
- direct solver 是否正常完成
- 是否触发孔隙度、体积分数、Biot 系数或 intrinsic/effective 输入检查

### 4.2 FIM vs Sequential 同源对比

比较：

- matrix pressure: `pressure_matrix_history.hdf5`
- fracture pressure: `pressure_fracture_history.hdf5`
- top displacement: `displacement_history.hdf5`

目标：

- 两个输入使用同一套同源 pressure 标定的初始压力和位移荷载
- 两个输入都使用 `crossStorageOffDiagScale="1.0"`
- 若两者差异异常，应先检查运行使用的 XML 副本、输出目录和绘图取样位置是否一致

### 4.3 任意 GEOS 案例 vs 解析解

比较：

- 初始压力 `p_m(0+)`、`p_f(0+)`
- Mandel-Cryer 早期峰值
- fracture 排水时间尺度
- matrix 中期平台
- matrix 晚期衰减
- top displacement

不同 GEOS 案例不要求逐点完全一致，但差异必须可解释为输入参数、耦合算法或数值离散差异，
而不是孔隙度体积分数处理、储量项或 Jacobian/残差实现不一致。判断结论应写入
`DPDP_Mandel_findings_current_status.md`，不要写入本 workflow。

当前阶段对 Sequential 的最低要求是：

- 量化记录裂缝排水时间尺度与手动 Fig. 5c 点的偏差
- `0.1s-10s`、`8s-30s`、以及 `tau=0.1-10` 范围内，中心点 matrix pressure 不应出现
  由外耦合提前停止造成的锯齿
- time step cut 为 0 或有明确原因
- 若 pressure 曲线与解析解不一致，不在本阶段强行调参补偿；该问题归入
  `DPDP_Mandel_findings_current_status.md` 的 FIM/Seq 精度问题

### 4.4 GEOS vs 手动解析点

压力主判据使用手动解析 CSV：

- matrix: `analitical_result/fig5c_primary_analitical.csv`
- fracture: `analitical_result/fig5c_secondary_analitical.csv`

GEOS 压力归一化：

- matrix: `p_m / 4.55e5`
- fracture: `p_f / 4.88e5`

时间归一化：

- `tau = t / t0`
- `t0` 使用 `script/dpdp_mandel_analytical.py` 中定义的同一值

推荐误差表采样点：

```text
1e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1, 3, 10, 30, 100, 300, 1000, 3000
```

误差表至少包含：

```text
tau,
manual_matrix,
fim_matrix,
seq_matrix,
fim_matrix_rel_error,
seq_matrix_rel_error,
manual_fracture,
fim_fracture,
seq_fracture,
fim_fracture_rel_error,
seq_fracture_rel_error
```

对 GEOS history 和手动解析点均采用插值到上述 `tau` 点的方式进行比较。若手动解析 CSV 在某些
`tau` 点外推风险较大，应在表中标注并避免把该点作为严格判据。

### 4.5 绘图脚本流程

通用绘图脚本：

```bash
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/script/plot_GEOS_vs_analytical.py
```

脚本不固定案例数量。默认压力基准为手动 Fig. 5c CSV：

- matrix: `analitical_result/fig5c_primary_analitical.csv`
- fracture: `analitical_result/fig5c_secondary_analitical.csv`

解析脚本曲线只作为可选诊断叠加，需要时加 `--show-script-analytical`。每个案例用
`--case "标签=输出目录"` 传入，输出目录必须包含：

- `pressure_matrix_history.hdf5`
- `pressure_fracture_history.hdf5`

显式指定案例的推荐命令：

```bash
python3 inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/script/plot_GEOS_vs_analytical.py \
  --case "FIM=<fim_sameSourcePressure_output_dir>" \
  --case "Seq=<seq_sameSourcePressure_output_dir>" \
  --show-script-analytical \
  --out <plot_output_dir>/geos_vs_fig5c_pressure.png
```

以后增加或减少案例时，只需增删 `--case` 参数。例如只比较两个 FIM，传两个 `--case`；
增加新的 Seq/FIM/诊断案例，则继续追加 `--case "标签=<output_dir>"`。

也可以让脚本递归搜索输出目录：

```bash
python3 inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/script/plot_GEOS_vs_analytical.py \
  --search-root inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs \
  --out inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/geos_vs_fig5c_pressure.png
```

自动搜索时，脚本会把所有包含 matrix/fracture pressure history 的目录都作为案例。标签优先从
输出目录中的 XML 副本推断；如果输出目录没有 XML，则使用目录名。因此推荐每次运行时把实际 XML
副本也放进输出目录。

脚本输出：

- PNG 图：`--out` 指定的文件
- 同名 PDF
- `<stem>_summary.csv`，包含每个案例的峰值、终值、以及相对手动 Fig. 5c CSV 的平均/最大绝对误差
- `<stem>_samples.csv`，包含固定 `tau` 采样点上的手动 CSV、每个 GEOS 案例和绝对误差

若把最终主对比图归档到 `analitical_result/`，文件名必须包含日期时间戳，例如：

```text
GEOS_pressure_sameSource_FIM_Seq_analytical_YYYYMMDD_HHMM.png
```

解析解脚本路径已经改为相对当前验证目录，不再依赖机器上的绝对路径。

### 4.6 当前 Seq 时间尺度检查记录

最近一次已完成的 Seq 时间尺度和平滑性检查：

- XML: `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml`
- 输出目录:
  `/tmp/dpdp_same_source_smooth_dt300_full_20260716_124345/seq`
- log:
  `/tmp/dpdp_same_source_smooth_dt300_full_20260716_124345/seq.log`
- 图:
  `analitical_result/GEOS_pressure_sameSource_FIM_Seq_analytical_20260716_1245_CN.png`
- 设置:
  `useIntrinsicInput=0`、fracture effective `K/G=4.514e8/3.108e8`、
  `maxSequentialPressureChange=100 Pa`、`10000 s` 之后 `forceDt=300 s`
- 运行结果:
  `1373` 个时间步，`0` 次 time step cut，总时间约 `29.8 s`

排水时间尺度：

```text
p_f/p0 = 0.5: manual 0.203984, Seq 0.421231
p_f/p0 = 0.1: manual 0.551803, Seq 1.080100
```

平滑性检查：

```text
actual 0.1s-10s: matrix pressure 一阶差分无符号反转
actual 8s-30s:   matrix pressure 一阶差分无符号反转
tau 0.1-10:      matrix pressure 一阶差分无符号反转
tau 800-1200:    matrix pressure 一阶差分无符号反转，最大斜率跳变约 0.0136
tau 9000-10500:  matrix pressure 一阶差分无符号反转，最大斜率跳变约 0.0100
```

Seq 求解策略历史对照：

```text
旧保守策略基线: 约 776 步，4974 次 nonlinear iteration，0 次切步，约 35.5 s
旧当前默认策略: 422 步，2006 次 nonlinear iteration，0 次切步，约 15.5 s
早期细时间步对照: 734 步，2996 次 nonlinear iteration，0 次切步，约 31.0 s
当前 late-time 平滑策略: 1373 步，3145 次 nonlinear iteration，0 次切步，约 29.8 s
```

早期细时间步对照保留 `maxSequentialPressureChange=100 Pa`，只把 `0.1-30 s` 的时间步恢复到旧细分段。
结果没有改善 fracture 排水 crossing：

```text
p_f/p0 = 0.5: current 0.421253, fine-early-dt 0.421902
p_f/p0 = 0.1: current 1.080309, fine-early-dt 1.109396
```

因此当前 Seq 不应再使用旧保守时间步作为默认验证策略；但为了图形平滑和晚期快速下降段精度，
`10000 s` 之后保留 `forceDt=300 s`。剩余的排水偏差应作为模型/分裂精度问题跟踪。

未解决问题的当前状态见 `DPDP_Mandel_findings_current_status.md`。

## 5. 通过标准

### 5.1 阶段性通过：可运行与可对比

当前阶段只要求：

- 同源 FIM、同源 Sequential 两个保留算例完整跑完，time step cut 为 0 或有明确说明
- 通用绘图脚本能把两个保留算例与解析解放在同一张图中，并输出误差汇总 CSV
- 曲线没有由外耦合提前停止导致的锯齿
- 两个保留输入均明确使用 `crossStorageOffDiagScale="1.0"`
- 不通过混用 intrinsic/effective K、Biot、孔隙度来“调曲线”

### 5.2 完整通过：FIM/Seq 物理验证

完整一轮验证通过需要同时满足：

- 同源 FIM、同源 Sequential 两个保留算例均完整跑完
- 两个保留输入均使用 `crossStorageOffDiagScale="1.0"`
- FIM 与手动解析 pressure 点在峰值、平台、排水时间尺度和晚期衰减趋势上合理一致
- Sequential 与手动解析 pressure 点趋势一致，若与 FIM 有差异，需要给出可解释原因
- 未出现孔隙度、体积分数、Biot 系数、裂缝参考孔隙度相对体积等输入概念混用

截至 2026-07-16，当前同源验证已达到 5.1 的可运行与可对比阶段；5.2 的解析一致性仍需继续跟踪。
当前 FIM 和 Seq 与解析解的偏差见 `DPDP_Mandel_findings_current_status.md`。

## 6. 更新算例时的检查清单

修改或新增验证算例后，提交前检查：

- 是否仍为 10x10 标准网格，除非明确做网格收敛性测试
- 是否没有覆盖手动解析 CSV
- 若 `script/dpdp_mandel_analytical.py` 或解析基准发生变化，是否同步重建 XML 中的顶部位移加载函数、
  初始位移表和相关说明
- 是否没有把运行结果或试错结论写入 `problem_description/`
- `mandel_input_tables/` 是否只保留当前 XML 实际引用的表；历史加载/诊断表若删除，是否记录到
  `DPDP_Mandel_trial_history_and_cleanup.md`
- `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` 是否保持为当前主 Seq 时间尺度验证 deck
- Seq effective 输入模式中 matrix/fracture K/G 和 Biot 是否均为 effective 参数，没有再混入
  fracture intrinsic 小刚度
- Seq 的 `maxSequentialPressureChange` 是否保持在 `1.0e2 Pa` 量级，避免外耦合单步提前停止
- FIM/Seq 同源输入是否使用同一套 pressure 标定表和位移荷载，且均为 `crossStorageOffDiagScale="1.0"`
- 绘图脚本是否用 `--case` 明确标注了每条曲线对应的输出目录
- Seq/FIM 耦合路径、输入模式或 cross-storage 处理若发生变化，是否同步更新
  `DPDP_Mandel_coupling_flow_Seq_vs_FIM.md`
- 物性、均匀化公式、storage 定义、孔隙度/渗透率体积分数约定或 XML 参数含义若发生变化，
  是否同步更新 `DPDP_Mandel_theory_input_reference.md`
- fracture porosity 是否表示裂缝介质自身孔隙度，而不是总体孔隙度
- 总体裂缝孔隙贡献是否通过 `fractureVolumeFraction * fracturePorosity` 理解和检查
- 裂缝渗透率是否按当前模型要求使用 bulk flux contribution
- 新 log、误差表、图像是否保存在新的验证输出目录

## 7. 文档更新要求

每轮验证或排查结束后：

- 新发现的问题、证据、是否解决、解决办法，写入 `DPDP_Mandel_findings_current_status.md`。
- 输入文件新增、删除、改名、保留原因，写入 `DPDP_Mandel_trial_history_and_cleanup.md`。
- 运行命令、绘图命令、输出目录规范变化，写入本文档。
- 论文原文、问题定义、边界条件来源变化，写入 `problem_description/`。
- 理论公式、物性参数定义、模型解释变化，写入 `DPDP_Mandel_theory_input_reference.md`。
- Seq/FIM 耦合执行流程、输入模式分支、cross-storage 处理或 relaxation 默认/用法变化，
  写入 `DPDP_Mandel_coupling_flow_Seq_vs_FIM.md`。

不要把文件清理记录写入 findings；不要把新的技术结论直接写入 workflow。
