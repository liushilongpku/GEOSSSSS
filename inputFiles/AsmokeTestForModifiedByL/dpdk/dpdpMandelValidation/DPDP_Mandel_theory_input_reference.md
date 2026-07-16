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
intrinsic 参数仍要保留给双孔隙 storage 公式使用。

---

## 3. Storage 与通量约定

### 3.1 双孔隙 storage matrix

Mehrabian-Abousleiman 恒应变 storage 写成：

```text
1/Mbar_ii = v_i * (1/M_i^intr + alpha_i^2 / K_i) - abar_i^2 / Kbar
1/Mbar_ij = -abar_i * abar_j / Kbar
1/M_i^intr = (alpha_i - phi_i) / Ks + phi_i * cf
```

在 GEOS 中，这部分由 `DualContinuumCrossFlow::assembleCouplingTerms` 使用 intrinsic 参数组装。
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

因此当前 fracture permeability 输入为：

```text
kappa_f = 0.03 * 4.935e-15 = 1.4805e-16 m2
```

`permeabilityComponents` 中应直接写入该 bulk flux contribution。初始化代码不再额外对
`permeabilityComponents` 乘 `v_i`，避免重复体积分数加权。

### 3.3 孔隙度与体积分数

XML 中 `defaultReferencePorosity` 表示 continuum-local intrinsic porosity：

```text
phi_m = 0.14
phi_f = 0.95
```

总体 REV 孔隙贡献通过体积分数理解：

```text
matrix pore contribution   = v_m * phi_m
fracture pore contribution = v_f * phi_f
```

不要把 `defaultReferencePorosity` 直接改成 `v_i * phi_i`。

---

## 4. 输入模式

`SinglePhaseDualContinuumPoromechanics` 支持两种 material input mode，由 `useIntrinsicInput` 控制。

### 4.1 `useIntrinsicInput="0"`：effective input

这是默认模式，也是当前保留两个主验证 XML 的模式。

在该模式下：

- matrix/fracture 本构中的 `K/G/Biot` 已经是 effective 值；
- `DualContinuumCrossFlow` 中仍需提供 intrinsic matrix/fracture 参数，用于 storage 公式；
- 求解器不会再次均匀化本构参数；
- 当前主 Seq `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` 不使用复合压力路径。

当前 effective input 的关键值：

```text
matrixSolid:   K=4.514e8, G=3.108e8
matrixPorosity defaultBiotCoefficient=0.382
fractureSolid: K=4.514e8, G=3.108e8   (当前 Seq)
fracturePorosity defaultBiotCoefficient=0.601
DualContinuumCrossFlow intrinsicMatrixBulkModulus=1.1e9
DualContinuumCrossFlow intrinsicMatrixBiot=0.9593
DualContinuumCrossFlow intrinsicFractureBulkModulus=2.25e7
DualContinuumCrossFlow intrinsicFractureBiot=0.9992
```

注：FIM effective deck 中 fracture 的 ElasticIsotropic K/G 不承担独立力学网格，但 fracture Biot、
fracture porosity 和 K_upf/K_pfu 映射仍需要保持与当前实现约定一致。具体 deck 以当前 XML 为准。

### 4.2 `useIntrinsicInput="1"`：intrinsic input

该模式要求本构直接填 intrinsic matrix/fracture 参数。当前保留主验证 XML 中没有 intrinsic deck，
但代码仍支持该输入模式。

FullyImplicit 路径：

- 初始化时调用 `computeEffectiveFromIntrinsic`；
- 用 Reuss 公式计算 `Kbar/Gbar/abar_i`；
- matrix mechanics 被覆写为 effective drained moduli；
- matrix Biot 被覆写为 `abar_m`；
- fracture mechanics Biot wrapper 写入 `abar_f`，用于 K_upf/K_pfu；
- fracture constitutive Biot 保持 intrinsic，避免污染 fracture porosity/storage；
- intrinsic 参数推送给 `DualContinuumCrossFlow`，供 storage 使用。

Sequential 路径：

- 初始化时不永久覆写本构；
- 运行时使用 composite-pressure mapping；
- 复合压力满足：

```text
alpha_m * p_eq = abar_m * p_m + abar_f * p_f
```

不要把 `useIntrinsicInput=1` 与 effective 本构参数混用，否则会触发输入模式检查或造成二次均匀化。

---

## 5. 当前保留输入文件

截至 2026-07-16，当前主验证只保留两个 XML：

| 输入文件 | 模式 | 用途 |
|---|---|---|
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml` | FIM effective input, `crossStorageOffDiagScale=1.0` | 当前同源全耦合法输入 |
| `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` | Sequential effective input, `crossStorageOffDiagScale=1.0` | 当前同源迭代耦合法输入 |

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
| `useIntrinsicInput` | `0` | `0`=effective input；`1`=intrinsic input |
| `autoInitializeStress` | `0` | 自动初始化有效应力；当前主 XML 不依赖它 |

### 6.2 `DualContinuumCrossFlow`

| 参数 | 默认 | 当前用途 |
|---|---:|---|
| `fractureVolumeFraction` | `-1` | storage 中的 REV fracture volume fraction；当前显式设为 `0.03` |
| `intrinsicMatrixBiot` | `-1` | matrix intrinsic Biot，用于 storage |
| `intrinsicMatrixBulkModulus` | `-1` | matrix intrinsic K，用于 storage |
| `intrinsicFractureBiot` | `-1` | fracture intrinsic Biot，用于 storage |
| `intrinsicFractureBulkModulus` | `-1` | fracture intrinsic K，用于 storage |
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
