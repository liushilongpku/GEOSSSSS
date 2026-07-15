# DPDP Mandel Seq/FIM 耦合执行流程

本文只说明当前 GEOS 自定义双连续介质孔隙力学求解器
`SinglePhaseDualContinuumPoromechanics`（类 `DualContinuumPoromechanicsSolverBase`）
中 Sequential 与 FullyImplicit 两条耦合路径的实现流程。

验证流程见 `DPDP_Mandel_validation_workflow.md`；当前问题状态、是否解决、最新曲线结论见
`DPDP_Mandel_findings_current_status.md`。本文不再记录阶段性试错结论，避免流程说明和验证结论混在一起。

代码位置：

- 主求解器：`src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp`
- 跨流/跨储存：`src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/SinglePhaseDualContinuum.cpp`

---

## 0. 公共数据结构与约定

- **mesh1**：基质（matrix）网格，同时承载力学；`SolidMechanicsLagrangianFEM` 只在 mesh1 上求解。
- **mesh2**：裂缝（fracture）网格，只承载流动，无独立力学自由度。
- 两套网格几何同位：mesh1 单元 k 与 mesh2 单元 k 占据同一物理位置，代表同一点处的两个连续介质。
- 自由度：位移 `u`（mesh1 节点）、基质压力 `p_m`（mesh1 单元）、裂缝压力 `p_f`（mesh2 单元）。
- 本验证的体积分数为 `v_m=0.97, v_f=0.03`；当前 N=2 Mandel 输入中
  `interporosityExchangeCoefficient=0`，matrix/fracture 之间无直接质量交换。
- 耦合方式由 XML `<NonlinearSolverParameters couplingType="...">` 选择：
  `Sequential` 或 `FullyImplicit`。

Sequential 固定应力所需字段不是在 `registerDataOnMesh` 处分流启用，而是在
`postInputInitialization()` 中、子求解器 `registerDataOnMesh` 之前启用：

- `solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate()`
- `flowSolver()->enableFixedStressPoromechanicsUpdate()`

该路径还要求：

- `subcycling="1"`
- `sequentialConvergenceCriterion="SolutionIncrements"`

---

## 1. Sequential（外迭代固定应力分裂）

Sequential 由 `CoupledSolver` 外层迭代驱动。每个时间步内交替求解力学子问题与双连续介质流动子问题，
并通过 `mapSolutionBetweenSolvers(domain, solverType)` 在子求解器之间传递状态。

### 1.1 单个外层迭代

```text
时间步开始
  implicitStepSetup 保存 _n 状态、压力状态等

Flow -> Mechanics 映射
  relaxSequentialFlowPressure
  flowSolver()->saveSequentialIterationState
  copyFracturePressureToMesh1
  若 useIntrinsicInput != 0:
    swapToCompositePressure
  否则:
    保持当前 p_m / p_f 表示，不使用复合压力
  updateBulkDensity

求解力学子问题
  SolidMechanics 只在 mesh1 上求解位移
  位移增量给出单元体应变增量 dEps_v

Mechanics -> Flow 映射
  若使用复合压力:
    restoreCompositePressure
  在 mesh1 上计算平均体应变增量
  更新 matrix porosity/permeability/fluid state
  将同一共享体应变增量映射给 fracture
  更新 fracture porosity/permeability/fluid state

求解 DualContinuumFlow
  matrix flow + fracture flow + cross-flow/cross-storage coupling
  当前 Sequential split 下 cross-storage offdiag 对当前外迭代滞后

检查外层收敛
  当前主 Seq 使用 SolutionIncrements，即检查压力外迭代增量
```

### 1.2 effective-input 与 intrinsic-input 的区别

当前主 Seq 输入文件 `DPDP_N2_dispdriven_seq_eff.xml` 使用：

```xml
useIntrinsicInput="0"
```

因此它是 **effective-input Sequential**。在该路径中：

- matrix/fracture 力学本构中的 `K/G/Biot` 已经是有效介质输入；
- 力学映射不使用 `p_eq` 复合压力；
- matrix 固定应力映射使用 `avgStress_m = K_m * dEps_v`；
- fracture 固定应力映射使用 `avgStress_f = K_f * dEps_v`；
- `DualContinuumCrossFlow` 中的 intrinsic 参数只用于双孔隙 storage 公式检查和组装。

只有 `useIntrinsicInput != 0` 时，Sequential 才使用复合压力技巧：

```text
alpha_m * p_eq = abar_m * p_m + abar_f * p_f
```

这个 intrinsic-input 路径会在力学前临时把 `p_m` 换成 `p_eq`，力学后再恢复真实 `p_m`。
不能把这个分支误写成当前主 Seq 的通用流程。

### 1.3 当前 Seq 稳定化要点

当前 Seq 的关键稳定化不是 pressure relaxation，也不是给物理参数加无意义缩放，而是在
fixed-stress Sequential split 下滞后 cross-storage offdiag 项：

- residual 中仍包含 offdiag storage 对上一外迭代压力的贡献；
- Jacobian 中只放入当前连续介质自己的 diagonal storage correction；
- 外迭代收敛后，offdiag 固定点在压力外迭代容差内与目标方程一致。

该处理针对的是分裂算法的收敛性。历史上 Seq full run 卡住/发散的主要原因就是 offdiag storage
隐式进入双流动子问题后，外迭代不再收缩。当前问题状态见 `DPDP_Mandel_findings_current_status.md`。

---

## 2. FullyImplicit（FIM，单体整装耦合）

FIM 不做外层子问题分裂。`DualContinuumPoromechanicsSolverBase::assembleSystem` 在一次装配中把
`u`、`p_m`、`p_f` 及其跨网格、跨连续介质块装入同一个雅可比系统，然后做单体 Newton 迭代。

FIM 不启用 Sequential fixed-stress 字段。若 `enableFimCrossStorage="1"`，双孔隙 cross-storage
correction 在 FIM 路径中开启。

### 2.1 一次 FIM assembleSystem 的主要步骤

```text
Step 0  mapFractureDataToMatrix
        把 p_f、alpha_f、fracture DOF 编号等映射到 mesh1 同位单元。

Step 1  matrix 单体孔隙力学核
        在 mesh1 上组装 u + p_m 的标准 SinglePhasePoromechanics 块。

Step 1.5 updateFracturePorosityFixedStress
        用 mesh1 位移增量计算共享体应变增量，更新 mesh2 fracture porosity、
        fracture mass 及 dMass/dP。

Step 2  K_upf：fracture pressure -> displacement residual/Jacobian
        组装 fracture pressure 对力学方程的贡献。

Step 2b K_pfu：displacement -> fracture mass Jacobian
        组装位移对 fracture mass balance 的一致 Jacobian。

Step 3  matrix face flux
        组装 mesh1 matrix 流动通量。

Step 4  fracture flow + dual-continuum cross coupling
        组装 mesh2 fracture 流动以及 matrix/fracture cross-storage correction。
```

`updateFracturePorosityFixedStress` 名称沿用 fixed-stress，但在 FIM 路径中它用于把同一体应变增量注入
fracture porosity，并由 `K_pfu` 保证单体 Newton 的一致性。单相路径中主要形式为：

```text
phi_f = phi_f,n + alpha_f * dEps_v + (alpha_f - phi_ref_f) / K_s * (p_f - p_f,n)
```

实际代码还会同步 fracture mass、`dMass/dP` 和相关应力诊断，因此上式只表示孔隙度更新的核心项。

### 2.2 FIM cross-storage correction

FIM 的 `DualContinuumCrossFlow::assembleCouplingTerms` 会将双孔隙储量补成 Mehrabian-Abousleiman
N=2 形式：

```text
1/Mbar_ii = v_i * (1/M_i^intr + alpha_i^2/K_i) - abar_i^2/Kbar
1/Mbar_ij = -abar_i * abar_j / Kbar * crossStorageOffDiagScale
```

其中 intrinsic 参数来自 `DualContinuumCrossFlow` 的 intrinsic 输入；effective Biot/模量用于
力学侧和当前有效介质输入。

`crossStorageOffDiagScale="1.0"` 是无经验 offdiag 缩放对照。`crossStorageOffDiagScale="0.911"`
是经验拟合输入，不应误写成无修正物理输入。关于 0.911 和当前解析基准不一致问题的定位，以
`DPDP_Mandel_findings_current_status.md` 为准。

FIM Newton 步长由 `fimNewtonRelaxation` 控制。代码默认值是 `0.5`；当前保留的两个 FIM XML 显式设置为
`0.7`。

---

## 3. 当前保留输入文件

截至 2026-07-15，本目录用于主验证的 XML 为：

| 输入文件 | 用途 |
|---|---|
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10.xml` | FIM 0.911，经验拟合输入 |
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10_noCorrection.xml` | FIM 1.0，无经验 offdiag 缩放对照 |
| `DPDP_N2_dispdriven_seq_eff.xml` | 当前 effective-input Sequential 输入 |

旧的 `correctLF`、stress-load、intrinsic/effective 试验 deck 不再作为当前主验证输入。历史试错文件的
清理记录见 `DPDP_Mandel_trial_history_and_cleanup.md`。

---

## 4. 对比速查

| 维度 | Sequential | FIM |
|---|---|---|
| 求解结构 | 外层 Picard/固定应力分裂，交替解力学和双流动 | 单体雅可比，`u/p_m/p_f` 同时 Newton |
| 力学网格 | 只在 mesh1 解位移 | 只在 mesh1 解位移 |
| fracture 力学作用 | 通过压力映射、共享体应变和 fracture porosity 更新进入流动 | 显式组装 K_upf/K_pfu，并更新 fracture porosity/mass |
| 当前主输入模式 | `useIntrinsicInput=0`，effective-input，不使用复合压力 | effective-input FIM |
| 复合压力 | 仅 intrinsic-input Sequential 使用 | 不使用 |
| cross-storage offdiag | Sequential split 下对外迭代滞后，保证固定点迭代收缩 | 直接进入单体 residual/Jacobian |
| relaxation | 当前验证不需要 pressure relaxation | XML 可显式设置 `fimNewtonRelaxation`；代码默认 0.5 |
| 当前代表 deck | `DPDP_N2_dispdriven_seq_eff.xml` | `DPDP_N2_dispdriven_fim_eff_direct_mesh10*.xml` |

---

## 5. 当前结果去向

本文不记录“哪条曲线最好”“过冲是否恢复”“误差多少”等验证结论。当前结论统一维护在：

- `DPDP_Mandel_findings_current_status.md`

绘图和运行流程统一维护在：

- `DPDP_Mandel_validation_workflow.md`

新增或修改耦合实现后，应同步检查本文是否仍与代码流程一致；新增验证结论则只写入 findings。
