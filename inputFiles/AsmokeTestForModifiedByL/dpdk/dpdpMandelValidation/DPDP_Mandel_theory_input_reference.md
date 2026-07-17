# DPDP Mandel 理论参数与输入约定

本文只维护双孔隙双渗透（DPDP, N=2）Mandel 验证所需的理论定义、参数约定、输入模式和实现接口。
运行流程见 `DPDP_Mandel_validation_workflow.md`；Seq/FIM 执行流程见
`DPDP_Mandel_coupling_flow_Seq_vs_FIM.md`；当前问题、证据、是否解决和最新验证结论见
`DPDP_Mandel_findings_current_status.md`。

本文不记录阶段性试错结论、曲线优劣排序、过冲是否恢复等结果性内容。

关键代码：

- `src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp`
- `src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/SinglePhaseDualContinuum.cpp`
- `src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumCrossFlow.{hpp,cpp}`
- 解析解脚本：`inputFiles/AsmokeTestForModifiedByL/dpdk/dpdpMandelValidation/script/dpdp_mandel_analytical.py`

---

## 1. 问题与物性

验证目标是 Mehrabian & Abousleiman (2014) 的 GOM shale N=2 Mandel 问题。计算域为
`0.03 m x 0.03 m x 0.03 m`，两个连续介质为 matrix 和 macro-fracture。

| 介质 | 排水体模量 K | 泊松比 | 剪切模量 G | 体积分数 v | 孔隙度 phi | 渗透率 k | intrinsic Biot |
|---|---:|---:|---:|---:|---:|---:|---:|
| matrix | `1.1e9 Pa` | `0.22` | `7.574e8 Pa` | `0.97` | `0.14` | `4.935e-21 m2` | `0.9593` |
| fracture | `2.25e7 Pa` | `0.22` | `1.549e7 Pa` | `0.03` | `0.95` | `4.935e-15 m2` | `0.9992` |

公共参数：

- grain bulk modulus: `Ks = 2.7e10 Pa`
- fluid bulk modulus: `Kf = 1.744e9 Pa`
- fluid compressibility: `cf = 5.734e-10 1/Pa`
- viscosity: `mu = 1.0e-3 Pa*s`
- direct matrix/fracture exchange: `interporosityExchangeCoefficient = 0`

当前验证输入采用位移驱动 Mandel 边界。具体运行、输出目录、绘图和判据由
`DPDP_Mandel_validation_workflow.md` 维护。

---

## 2. 有效介质参数

DPDP Mandel 解析解使用共享体应变的有效介质。对于 N=2：

```text
Kbar = 1 / (v_m / K_m + v_f / K_f)
Gbar = 1 / (v_m / G_m + v_f / G_f)
abar_i = Kbar * v_i * alpha_i / K_i
```

当前物性给出：

```text
Kbar  ~= 4.514e8 Pa
Gbar  ~= 3.108e8 Pa
abar_m ~= 0.382
abar_f ~= 0.601
```

这些是 effective-input deck 中写入 matrix/fracture 力学本构和 BiotPorosity 的有效参数。
双孔隙 storage 可以用 direct effective storage 输入，避免在 effective-input deck 中继续混入
`intrinsic*` 参数名。

---

## 3. Storage 与通量约定

### 3.1 双孔隙 storage matrix

Mehrabian-Abousleiman 恒应变 storage 写成：

```text
1/Mbar_ii = v_i * (1/M_i^intr + alpha_i^2 / K_i) - abar_i^2 / Kbar
1/Mbar_ij = -abar_i * abar_j / Kbar
1/M_i^intr = (alpha_i - phi_i) / Ks + phi_i * cf
```

在 GEOS 中，这部分由 `DualContinuumCrossFlow::assembleCouplingTerms` 组装。当前支持两种一致输入方式：

- direct effective storage：直接给出 `effectiveMatrixStorage`、`effectiveFractureStorage`、
  `effectiveCrossStorage`；
- intrinsic reconstruction：`useIntrinsicInput="1"` 时从 matrix/fracture 材料本构读取本征
  `K/G/Biot/phi/Ks`，代码先均匀化 mechanics/Biot，再按上式重构 storage。

`DualContinuumCrossFlow` 不再暴露 `intrinsicMatrixBiot`、`intrinsicMatrixBulkModulus`、
`intrinsicFractureBiot`、`intrinsicFractureBulkModulus` XML 属性；这些本征值必须来自材料本构。
不要把 effective K/Biot 写进本征输入模式。

`crossStorageOffDiagScale` 只作用于 offdiag 项：

```text
1/Mbar_ij -> 1/Mbar_ij * crossStorageOffDiagScale
```

当前约定：

- `crossStorageOffDiagScale="1.0"`：无经验 offdiag 缩放，对照输入。
- `crossStorageOffDiagScale="0.911"`：历史经验拟合输入，已删除，不再作为当前验证输入。
- `0.911` 不应写成“无修正物理输入”；其当前定位和解析基准不一致问题以
  `DPDP_Mandel_findings_current_status.md` 为准。

### 3.2 渗透率与体积分数

双连续介质对 bulk Darcy flux 的贡献按体积分数加权：

```text
kappa_i = v_i * k_i
```

渗透率输入不由 `useIntrinsicInput` 切换，XML 始终写 continuum-local intrinsic `k_i`，代码将
permeability 乘 `v_i`。因此当前主 XML 中的 permeability 输入为：

```text
k_m = 4.935e-21 m2  -> REV contribution 0.97 * k_m = 4.78695e-21 m2
k_f = 4.935e-15 m2  -> REV contribution 0.03 * k_f = 1.4805e-16 m2
```

### 3.3 孔隙度与体积分数

孔隙度输入也不由 `useIntrinsicInput` 切换。XML 中 `defaultReferencePorosity` 始终表示
continuum-local intrinsic porosity：

```text
phi_m = 0.14
phi_f = 0.95
```

代码会乘体积分数得到总体 REV 孔隙贡献：

```text
matrix pore contribution   = v_m * phi_m
fracture pore contribution = v_f * phi_f
current REV contribution: phi_m_REV=0.1358, phi_f_REV=0.0285
```

---

## 4. 输入模式

`SinglePhaseDualContinuumPoromechanics` 支持两种 mechanics/storage input mode，由
`useIntrinsicInput` 控制。孔隙度和渗透率不参与该 flag 切换，始终是本征输入。

### 4.1 `useIntrinsicInput="0"`：effective mechanics/storage input

这是当前保留的 effective 验证 XML 显式使用的模式，但不再是代码默认模式。代码默认
`useIntrinsicInput="1"`，即 intrinsic mechanics/storage input。

在该模式下：

- matrix/fracture 本构中的 `K/G/Biot` 已经是 effective 值；
- matrix/fracture 的 `defaultReferencePorosity` 和 `permeabilityComponents` 仍是 intrinsic 值，
  代码会乘 `v_m/v_f`；
- 当前主 XML 用 `effectiveMatrixStorage/effectiveFractureStorage/effectiveCrossStorage` 直接输入
  双孔隙 storage matrix；
- 不允许提供旧 CrossFlow `intrinsic*` storage reconstruction 属性；schema 中已删除这些属性；
- 求解器不会再次均匀化本构参数；
- 当前主 Seq `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` 不使用复合压力路径。

当前 effective input 的关键值：

```text
matrixSolid:   K=4.514e8, G=3.108e8
matrixPorosity defaultReferencePorosity=0.14
matrixPorosity defaultBiotCoefficient=0.382
matrixPerm permeabilityComponents=4.935e-21
fractureSolid: K=4.514e8, G=3.108e8
fracturePorosity defaultReferencePorosity=0.95
fracturePorosity defaultBiotCoefficient=0.601
fracturePerm permeabilityComponents=4.935e-15
DualContinuumCrossFlow effectiveMatrixStorage=5.95755344271165954e-10
DualContinuumCrossFlow effectiveFractureStorage=5.46327412521659742e-10
DualContinuumCrossFlow effectiveCrossStorage=-5.08769676169630689e-10
```

### 4.2 `useIntrinsicInput="1"`：intrinsic mechanics/storage input

这是代码默认模式。该模式要求本构直接填 intrinsic matrix/fracture 参数，包括 intrinsic
`defaultReferencePorosity` 和 intrinsic permeability。不允许提供 direct effective storage；代码会对
孔隙度和渗透率乘 `v_m/v_f`，并由本征参数重构 storage。当前保留 FIM/Seq intrinsic
deck 用于检查该路径与 effective 输入路径的一致性。

FullyImplicit 与 Sequential 路径：

- 初始化时调用 `computeEffectiveFromIntrinsic`；
- 用 Reuss 公式计算 `Kbar/Gbar/abar_i`；
- matrix/fracture mechanics 被覆写为 effective drained moduli；
- matrix/fracture Biot 被覆写为 `abar_m/abar_f`；
- fracture mechanics Biot wrapper 写入 `abar_f`，用于 K_upf/K_pfu；
- intrinsic 参数推送给 `DualContinuumCrossFlow`，供 storage 使用。

历史 Sequential runtime composite-pressure mapping 已禁用；受控对比显示，只临时替换 matrix
pressure 而不同时均匀化 fracture 本构，会造成 intrinsic Seq 与 effective Seq 不一致。

不要把 `useIntrinsicInput=1` 与 effective 本构参数混用，否则会触发输入模式检查或造成二次均匀化。

---

## 5. 当前保留输入文件

截至 2026-07-17，当前验证保留四个 XML：

| 输入文件 | 模式 | 用途 |
|---|---|---|
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml` | FIM effective input, `crossStorageOffDiagScale=1.0` | 当前同源全耦合法输入 |
| `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` | Sequential effective input, `crossStorageOffDiagScale=1.0` | 当前同源迭代耦合法输入 |
| `DPDP_N2_dispdriven_fim_intrinsic_sameSourcePressure.xml` | FIM intrinsic input, `crossStorageOffDiagScale=1.0` | 本征输入一致性检查 |
| `DPDP_N2_dispdriven_seq_intrinsic_sameSourcePressure.xml` | Sequential intrinsic input, `crossStorageOffDiagScale=1.0` | 本征输入一致性检查 |

历史 intrinsic、stress-load、correctLF、kappa 扫描等 deck 的保留/删除原因见
`DPDP_Mandel_trial_history_and_cleanup.md`。

---

## 6. XML 参数速查

### 6.1 `SinglePhaseDualContinuumPoromechanics`

| 参数 | 默认 | 当前用途 |
|---|---:|---|
| `fractureVolumeFraction` | `-1` | REV fracture volume fraction；当前显式设为 `0.03` |
| `fimNewtonRelaxation` | `0.5` | FIM Newton 欠松弛；当前 FIM XML 显式设为 `0.7` |
| `sequentialPressureRelaxation` | `1.0` | Seq pressure under-relaxation；当前验证不需要 relaxation |
| `enableFractureMechanicsCoupling` | `1` | FIM 中装配 K_upf/K_pfu |
| `enableFracturePorosityStrainCoupling` | `1` | FIM 中启用 strain -> fracture porosity 及 K_pfu |
| `enableFimCrossStorage` | `1` | FIM 中启用双孔隙 cross-storage correction |
| `useIntrinsicInput` | `1` | `1`=intrinsic mechanics/storage 并由代码均匀化与重构 storage；`0`=effective mechanics/storage 且必须提供 direct effective storage；porosity/permeability 始终为本征输入 |
| `autoInitializeStress` | `0` | 自动初始化有效应力；当前主 XML 不依赖它 |

### 6.2 `DualContinuumCrossFlow`

| 参数 | 默认 | 当前用途 |
|---|---:|---|
| `fractureVolumeFraction` | `-1` | storage 中的 REV fracture volume fraction；当前显式设为 `0.03` |
| `effectiveMatrixStorage` | `0` | direct effective matrix storage `Sbar_mm`；当前主 XML 使用 |
| `effectiveFractureStorage` | `0` | direct effective fracture storage `Sbar_ff`；当前主 XML 使用 |
| `effectiveCrossStorage` | `0` | direct effective offdiag storage `Sbar_mf=Sbar_fm`；当前主 XML 使用 |
| `crossStorageOffDiagScale` | `1.0` | offdiag storage scale；`1.0` 为无经验缩放对照 |
| `interporosityExchangeCoefficient` | `0` | N=2 Mandel 当前取 `0` |

修改参数说明字符串后，需要重新生成 schema：

```bash
make -C build geosx_generate_schema
```

---

## 7. 文档边界

本文只维护技术定义和输入约定：

- 物性、体积分数、孔隙度、渗透率加权
- effective/intrinsic 参数关系
- XML 参数含义
- 解析解脚本路径和主基准数据说明

以下内容不写入本文：

- 最新运行是否通过
- FIM/Seq 哪个拟合更好
- 过冲、排水、晚期衰减等问题是否解决
- 删除/保留历史 XML 的文件操作记录

这些内容分别维护在：

- 问题与结论：`DPDP_Mandel_findings_current_status.md`
- 运行与绘图流程：`DPDP_Mandel_validation_workflow.md`
- 文件清理记录：`DPDP_Mandel_trial_history_and_cleanup.md`
- Seq/FIM 执行顺序：`DPDP_Mandel_coupling_flow_Seq_vs_FIM.md`
