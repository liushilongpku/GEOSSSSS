# DPDP Mandel FIM/Seq 验证流程

本文档用于在双重介质流固耦合程序改动后，重新验证本目录中的 FIM 与 Sequential
Mandel 算例。目标是让后续更新验证算例时有固定流程，避免混用历史试验文件、覆盖旧结果，
或把阶段性调通结果误认为最终物理验证通过。

当前阶段的工作顺序：

1. 先把排水时间尺度和数值曲线平滑性大致对上。
2. 再处理 Mandel-Cryer 过冲、加载方式和基质早期响应。
3. 最后做 FIM effective、FIM intrinsic、Sequential 三者的完整一致性验证。

## 1. 验证基准

### 1.1 解析解基准

`analitical_result/fig5c_primary_analitical.csv` 是手动取得的 Fig. 5c primary/matrix
压力解析解点。

`analitical_result/fig5c_secondary_analitical.csv` 是手动取得的 Fig. 5c secondary/fracture
压力解析解点。

这两个 CSV 是本验证流程的压力主基准，不应当被脚本重新生成或覆盖。脚本解析解
`script/dpdp_mandel_analytical.py` 可作为辅助检查，但最终误差表应以这两个手动 CSV 为准。

`analitical_result/fig5d_analitical.csv` 可作为应力辅助诊断。当前位移驱动算例中，应力需要由
有效应力和 Biot 压力项重构，不作为主要通过/失败判据。

### 1.2 网格

使用 10x10 网格作为标准验证网格。当前已经认为 10x10 满足该 Mandel 验证精度要求，无需使用
20x20 或 40x40 作为常规验证。

## 2. 推荐算例

### 2.1 FIM 有效输入主算例

基准文件：

`DPDP_N2_dispdriven_fim_eff_direct_mesh10.xml`

验证前必须确认：

- `couplingType="FullyImplicit"`
- `enableFimCrossStorage="1"`
- `crossStorageOffDiagScale="1.0"`
- mesh 为 `nx="{ 10 }"`、`nz="{ 10 }"`
- 裂缝体积分数 `fractureVolumeFraction="0.03"`
- 裂缝渗透率采用体积加权后的 bulk flux contribution，即 `kappa_f = v_f * k_f`
- `defaultReferencePorosity` 输入 intrinsic 孔隙度；程序只将 `v_i` 乘入 storage/reference porosity。
  `permeabilityComponents` 不会再由 dual-continuum 初始化乘 `v_i`，因此本验证中必须直接输入 `kappa_i`
- 不使用历史调参值 `crossStorageOffDiagScale="0.911"`；该值来自旧储量项错误阶段的补偿，不作为当前验证设置

### 2.2 FIM intrinsic 输入一致性算例

基准文件：

`DPDP_N2_dispdriven_fim_eff_direct_mesh10_INTRINSIC.xml`

验证前必须确认：

- `couplingType="FullyImplicit"`
- `useIntrinsicInput="1"`
- `enableFimCrossStorage="1"`
- `crossStorageOffDiagScale="1.0"`
- intrinsic 材料输入为 `K_m=1.1e9`、`alpha_m=0.9593`、`K_f=2.25e7`、`alpha_f=0.9992`
- 自动均匀化后的结果应与 FIM 有效输入主算例基本一致

### 2.3 Sequential 对照算例

当前优先使用：

`DPDP_N2_dispdriven_seq_eff.xml`

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
- `maxSequentialPressureChange="1.0e1"`。不要恢复到 `1.0e3`；`1.0e3` 会导致部分时间步
  外耦合只迭代一次、下一步再补偿，在 `tau≈0.8-1` 附近形成压力锯齿
- 裂缝渗透率采用 `kappa_f = v_f * k_f`
- 孔隙度输入为 intrinsic 孔隙度，裂缝总体孔隙贡献按
  `fractureVolumeFraction * fracturePorosity` 理解和检查

历史文件 `DPDP_N2_dispdriven_correctLF_kappa.xml` 可作为加载函数/渗透率处理的对照，但当前
Seq 时间尺度验证以 `DPDP_N2_dispdriven_seq_eff.xml` 为准。

## 3. 运行规范

### 3.1 输出目录

每次验证必须使用新的独立输出目录，避免覆盖历史结果。推荐格式：

```bash
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<case>_<date>_<tag>/
```

每个算例单独建子目录，例如：

```bash
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/fim_eff_mesh10_YYYYMMDD_01/
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/fim_intrinsic_mesh10_YYYYMMDD_01/
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/seq_full_effective_fracKeff_tol10_single_YYYYMMDD_01/
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
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/DPDP_N2_dispdriven_seq_eff.xml \
  -o inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<seq_output_dir> \
  > inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/runs/<seq_output_dir>.log 2>&1
```

MPI 命令模板：

```bash
mpirun -np 4 build-ubuntu-lsl-release/bin/geosx -x 2 -y 2 -z 1 \
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

### 4.2 FIM effective vs FIM intrinsic

比较：

- matrix pressure: `pressure_matrix_history.hdf5`
- fracture pressure: `pressure_fracture_history.hdf5`
- top displacement: `displacement_history.hdf5`

目标：

- 两者在相同时间点上的曲线基本重合
- 若有差异，应优先检查 `useIntrinsicInput` 自动均匀化是否与手工 effective 输入一致
- 特别检查等效 Biot 系数、等效骨架模量、裂缝参考孔隙度与体积分数是否被重复加权或漏加权

### 4.3 FIM vs Sequential

比较：

- 初始压力 `p_m(0+)`、`p_f(0+)`
- Mandel-Cryer 早期峰值
- fracture 排水时间尺度
- matrix 中期平台
- matrix 晚期衰减
- top displacement

FIM 与 Sequential 不要求逐点完全一致，但差异必须可解释为耦合算法差异，而不是输入参数、孔隙度体积分数处理、储量项或 Jacobian/残差实现不一致。

当前阶段对 Sequential 的最低要求是：

- 裂缝排水时间尺度与手动 Fig. 5c 点基本一致
- `0.1s-10s`、`8s-30s`、以及 `tau=0.1-10` 范围内，中心点 matrix pressure 不应出现
  由外耦合提前停止造成的锯齿
- time step cut 为 0 或有明确原因
- 若 matrix 早期响应和解析解不一致，不在本阶段强行调参补偿；该问题归入后续
  Mandel-Cryer 过冲/加载方式验证

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

### 4.5 当前 Seq 时间尺度检查记录

最近一次已完成的 Seq 时间尺度检查：

- XML: `DPDP_N2_dispdriven_seq_eff.xml`
- 输出目录:
  `runs/seq_full_effective_fracKeff_tol10_single_20260714_01/`
- log:
  `runs/seq_full_effective_fracKeff_tol10_single_20260714_01.log`
- 图:
  `runs/seq_full_effective_fracKeff_tol10_single_20260714_01/seq_drainage_time_comparison_fracKeff_tol10_20260714.png`
- 设置:
  `useIntrinsicInput=0`、fracture effective `K/G=4.514e8/3.108e8`、
  `maxSequentialPressureChange=10 Pa`
- 运行结果:
  `776` 个时间步，`0` 次 time step cut，总时间约 `16.85 s`

排水时间尺度：

```text
p_f/p0 = 0.5: manual 0.203984, Seq 0.224645
p_f/p0 = 0.1: manual 0.551803, Seq 0.634275
```

平滑性检查：

```text
actual 0.1s-10s: matrix pressure 一阶差分无符号反转
actual 8s-30s:   matrix pressure 一阶差分无符号反转
tau 0.1-10:      matrix pressure 一阶差分无符号反转
```

未解决问题：

- matrix 早期响应和手动解析解仍不一致
- Mandel-Cryer 过冲尚未恢复
- 当前位移驱动加载方式与常载荷解析解之间的对应关系需要后续单独验证

## 5. 通过标准

### 5.1 阶段性通过：排水时间尺度

当前阶段只要求：

- Seq 算例完整跑完，time step cut 为 0 或有明确说明
- fracture 排水时间尺度与手动解析点在同一量级并接近
- 曲线没有由外耦合提前停止导致的锯齿
- 不使用 `crossStorageOffDiagScale="0.911"` 等历史补偿参数
- 不通过混用 intrinsic/effective K、Biot、孔隙度来“调曲线”

### 5.2 完整通过：FIM/Seq 物理验证

完整一轮验证通过需要同时满足：

- FIM effective、FIM intrinsic、Sequential 三个主算例均完整跑完
- FIM 使用 `crossStorageOffDiagScale="1.0"`，不依赖历史补偿值 `0.911`
- FIM effective 与 FIM intrinsic 的压力和位移曲线基本一致
- FIM 与手动解析 pressure 点在峰值、平台、排水时间尺度和晚期衰减趋势上合理一致
- Sequential 与手动解析 pressure 点趋势一致，若与 FIM 有差异，需要给出可解释原因
- 未出现孔隙度、体积分数、Biot 系数、裂缝参考孔隙度相对体积等输入概念混用

## 6. 更新算例时的检查清单

修改或新增验证算例后，提交前检查：

- 是否仍为 10x10 标准网格，除非明确做网格收敛性测试
- 是否没有覆盖手动解析 CSV
- `DPDP_N2_dispdriven_seq_eff.xml` 是否保持为当前主 Seq 时间尺度验证 deck
- Seq effective 输入模式中 matrix/fracture K/G 和 Biot 是否均为 effective 参数，没有再混入
  fracture intrinsic 小刚度
- Seq 的 `maxSequentialPressureChange` 是否保持在 `1.0e1 Pa` 量级，避免外耦合单步提前停止
- FIM 主算例是否使用 `crossStorageOffDiagScale="1.0"`
- intrinsic 输入算例是否开启 `useIntrinsicInput="1"`
- fracture porosity 是否表示裂缝介质自身孔隙度，而不是总体孔隙度
- 总体裂缝孔隙贡献是否通过 `fractureVolumeFraction * fracturePorosity` 理解和检查
- 裂缝渗透率是否按当前模型要求使用 bulk flux contribution
- 新 log、误差表、图像是否保存在新的验证输出目录
