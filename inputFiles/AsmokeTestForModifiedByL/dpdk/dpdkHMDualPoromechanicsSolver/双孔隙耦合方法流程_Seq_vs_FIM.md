# 双孔隙双渗透（DPDP, N=2）孔隙力学耦合方法详解：Sequential vs FIM

本文整理 GEOS 自定义双连续介质孔隙力学求解器
`SinglePhaseDualContinuumPoromechanics`（类 `DualContinuumPoromechanicsSolverBase`）
两种耦合策略的**详细执行流程**，对应 Mehrabian & Abousleiman (2014) 的 N=2 Mandel 问题。

代码位置：
- 主求解器：`src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp`
- 跨流/跨储存：`src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/SinglePhaseDualContinuum.cpp`

---

## 0. 公共数据结构与约定

- **mesh1**：基质（matrix）网格，**同时承载力学**（`SolidMechanicsLagrangianFEM` 只在 mesh1 上求解）。
- **mesh2**：裂缝（fracture）网格，只有流动，无独立力学。
- 两套网格**几何同位（co-located）**：mesh1 单元 k 与 mesh2 单元 k 占据同一物理位置，代表同一点处的两个连续介质。
- 自由度：位移 `u`（mesh1 节点）、基质压力 `p_m`（mesh1 单元）、裂缝压力 `p_f`（mesh2 单元）。
- 体积分数 `v_m=0.97, v_f=0.03`；耦合通过**共享体应变**实现（本算例 `interporosityExchangeCoefficient=0`，无直接质量交换）。
- 耦合方式由 XML `<NonlinearSolverParameters couplingType="...">` 选择：`Sequential` 或 `FullyImplicit`。

两条路径在 `registerDataOnMesh` 处就分流：仅当 `couplingType==Sequential` 时调用
`enableFixedStressPoromechanicsUpdate()`（注册 `pressure_k` 等固定应力所需场）。

---

## 1. Sequential（迭代/固定应力耦合）

外层由基类 `CoupledSolver` 的顺序迭代驱动；每个时间步内做**外层 Picard 迭代**，交替求解力学子问题与流动子问题，直到外层收敛。子求解器之间的状态传递全部由本类重写的
`mapSolutionBetweenSolvers(domain, solverType)` 完成。

### 1.1 单个外层迭代的执行顺序

```
┌─ 时间步开始（implicitStepSetup：保存 _n 状态、p_eq_n 等）
│
├─ [子问题 A] 求解力学（SolidMechanics，仅 mesh1）
│     ├─ 进入前 mapSolutionBetweenSolvers(Flow) 已为力学准备好应力/压力：
│     │     · copyFracturePressureToMesh1：把 p_f 拷到 mesh1 单元场
│     │     · swapToCompositePressure：令力学看到“复合压力” p_eq，
│     │       满足 α_m·p_eq = ᾱ_m·p_m + ᾱ_f·p_f（用有效 Biot 合成单一等效压力）
│     │     · updateBulkDensity
│     └─ 力学解出位移增量 du → 单元体应变增量 dEps_v、平均平均总应力增量
│
├─ mapSolutionBetweenSolvers(SolidMechanics)：力学→流动 传递
│     ├─ restoreCompositePressure：恢复真实 p_m（撤销 swap）
│     ├─ 在 mesh1 上对每单元做平均平均总应力增量
│     │   AverageOverQuadraturePoints（mesh2 无 FE，不做）
│     └─ 固定应力孔隙度更新（clean dual-continuum fixed-stress）：
│         力学解的是“复合介质”，有效排水模量 K_eff、复合压力 p_eq：
│           dSigma_v = K_eff·dEps_v − α_m·dp_eq
│         共享体应变   dEps_v = (dSigma_v + α_m·dp_eq)/K_eff
│         给每个连续介质 i 注入其有效-Biot 应变耦合 ᾱ_i·dEps_v，
│         做法：令 avgStress_i = v_i·K_eff·dEps_v，则 BiotPorosity 的
│         α_i·avgStress_i/K_i 恰好 = ᾱ_i·dEps_v（矩阵用自身 K_m、α_m，
│         与裂缝对称，无 K 来回交换）。
│
├─ [子问题 B] 求解流动（DualContinuumFlow = 基质流 + 裂缝流 + 跨流耦合）
│     ├─ 基质流 SinglePhaseFVM（mesh1）：储存项含上一步固定应力孔隙度
│     │   d φ/dp = (α−φ)/Ks + φ·c_f + α²/K_fs（固定应力稳定项）
│     ├─ 裂缝流 SinglePhaseFVM（mesh2）
│     └─ DualContinuumCrossFlow::assembleCouplingTerms：跨连续介质质量交换
│         （本算例 Γ=0，主要是 Sequential 路径下的多孔隙储存修正）
│
└─ 检查外层收敛（pressure_k 的变化 < 容差）；未收敛则回到子问题 A
   收敛后 implicitStepComplete（推进 _n 状态）
```

### 1.2 Sequential 的核心思想与已知问题

- **固定应力分裂（fixed-stress split）**：流动子问题里把力学的反馈用一个**滞后**的
  稳定项 `α²/K_fs` 近似进储存（`computePorosityFixedStress`），力学子问题里用滞后压力。
  两个子问题各自良态、可分别用各自的线性求解器。
- **复合压力技巧**：力学只在 mesh1 上、只解一套位移，却要同时感受 p_m 与 p_f 的总应力贡献，
  因此用 `α_m·p_eq = ᾱ_m·p_m + ᾱ_f·p_f` 合成单一等效压力 p_eq 喂给力学，解后再恢复。
- **已知问题（软裂缝下排水偏慢）**：固定应力稳定项 `α_f²/K_f`（软裂缝 K_f=2.25e7）
  把裂缝储存放大约 80×；若外层 subcycling 不足（单趟），这个放大未被充分抵消，
  导致排水过程偏慢。详见 memory `dpdp-mandel-oscillation`。

---

## 2. FullyImplicit（FIM，单体整装耦合）

不做子问题分裂；本类**重写 `assembleSystem`**，在一次装配里把 u、p_m、p_f
所有块（含跨网格、跨连续介质）装进同一个雅可比，然后做单体 Newton 迭代。
`mapSolutionBetweenSolvers` 的 Sequential 逻辑在 FIM 下不走（固定应力不启用）。

### 2.1 一次 assembleSystem 的装配步骤

```
Step 0  mapFractureDataToMatrix
        把 p_f、α_f、裂缝 DOF 编号从 mesh2 映射到 mesh1 的同位单元，
        供后续跨网格块按行/列正确寻址。

Step 1  单体基质核（u + p_m）
        在 mesh1 上跑成熟的 SinglePhasePoromechanics 核，一次积分循环产出
        K_uu、K_upm、K_pmu、K_pmpm（同一本构状态 → Schur 补自洽，
        对超低渗也稳定）。

Step 1.5 updateFracturePorosityFixedStress（名为 fixed-stress，实为真实 Biot 更新）
        用 mesh1 的位移增量算体应变，更新裂缝孔隙度与裂缝质量及其对 p_f 的导数：
          φ_f = φ_f,n + α_f·dEps_v + (α_f−φ_f)/Ks·Δp_f
        （储存用真实 (α_f−φ)/Ks，无 α²/K 放大；应变耦合显式带入）。

Step 2  K_upf：裂缝压力 → 位移残差
        +α_f·p_f·∇N（总应力约定）。跨网格 u↔p_f 稀疏性由
        DofManager::addCouplingDualContinuumMechanics 提供，条目不再被丢弃。

Step 2b K_pfu：位移 → 裂缝质量雅可比
        d(φ_f·ρ·V)/dU = ρ·α_f·∇N，使单体 Newton 对 u↔p_f 一致。

Step 3  基质面通量 assembleFluxTerms（mesh1 TPFA）。

Step 4  裂缝流（SinglePhaseFVM/secondarySolver，mesh2）+ 跨流。

Step 4b DualContinuumCrossFlow::assembleCouplingTerms 的多孔隙储存修正：
        把每单元储存补成 Mehrabian eq.A25 的 1/M̄ 矩阵：
          对角  1/M̄_ii = v_i(1/M_i^intr + α_i²/K_i) − ᾱ_i²/K̄
          非对角 1/M̄_ij = −ᾱ_i·ᾱ_j/K̄ × crossStorageOffDiagScale
```

装配后求解线性系统、用 `scalingForSystemSolution` 施加 **Newton 欠松弛**
（`fimNewtonRelaxation`，默认 0.7），抑制软耦合下的 period-2（λ≈−1）振荡，
迭代至单体 Newton 收敛。

### 2.2 FIM 的核心思想与已知问题

- **整装雅可比**：u、p_m、p_f 全在一个系统里，Schur 补自洽 → 收敛稳健，无固定应力滞后。
- **有效介质参数**：mesh1 用排水有效模量 K̄=4.514e8、Ḡ=3.108e8、有效 Biot ᾱ=[0.382,0.601]；
  储存按 eq.A25 装配。位移 BC 用解析 u_z(t)（位移驱动）。
- **已知问题（非对角近抵消，需 0.911）**：有效扩散储存 `Sbar = invM + ᾱ⊗cm`。
  FIM 把全量 `invM` 非对角 `−ᾱ_m·ᾱ_f/K̄` 放进储存，靠力学 Schur 补 `+ᾱ_m·ᾱ_f/(K̄+4Ḡ/3)`。
  这是两个 ~5e-10 大数的**近抵消**；实测离散 Q1 力学比连续 oedometric **刚 21%**，
  抵消不足 → 平台偏低。`crossStorageOffDiagScale=0.911`（实测、网格无关常数）修正之，
  平台对解析吻合到 <0.25%。详见 memory `dpdp-fim-crossstorage-offdiag` 与
  `SinglePhaseDualContinuum.cpp` 内注释。残留：Mandel-Cryer 过冲峰偏低约 0.10（独立瞬态问题）。

---

## 3. 对比速查

| 维度 | Sequential | FIM |
|---|---|---|
| 求解结构 | 外层 Picard 交替解 力学/流动 两个子系统 | 单次装配的整装雅可比 + 单体 Newton |
| 力学反馈 | 固定应力稳定项 α²/K_fs（滞后近似） | K_upf/K_pfu 真实雅可比块 |
| 子求解器 | 力学、流动各自独立线性求解 | 一个耦合线性系统 |
| 压力传递 | 复合压力 p_eq swap/restore | 直接用 p_m、p_f（Step 0 映射） |
| 孔隙度更新 | computePorosityFixedStress（含 α²/K） | 真实 Biot（Step 1 核 + Step 1.5） |
| 稳定/收敛手段 | 外层 subcycling | scalingForSystemSolution 欠松弛 0.7 |
| 主要短板 | 软裂缝 α²/K 放大 → 排水偏慢 | 非对角近抵消 → 需 0.911 标定；过冲略低 |
| 代表算例 | `DPDP_N2_dispdriven_correctLF.xml` | `DPDP_N2_dispdriven_fim_eff*.xml` |

---

## 4. 与解析解的吻合现状（位移驱动 Mandel, N=2）

- **FIM（fim_eff + scale=0.911）**：与验证过的解析解全程吻合 —— 平台 ~0.885 精确、排水跟随；
  仅过冲峰略低（~1.05 vs 1.08）。见 `analitical_result/GEOS_vs_analytical_vs_digitized_full.png`。
- **Sequential（cross-storage 版）**：矩阵平台 ~0.93、几乎无过冲（早期 cross-storage 里程碑，
  见 `analitical_result/GEOS_crossStorage_vs_analytical_vs_digitized.png`，为 6/2 的 Sequential 结果）。

> 注：本仓库解析脚本 `dpdkHMValidation/dpdp_mandel_analytical.py` 与论文数字化曲线在
> **晚期排水时间尺度**上存在 ~10× 差异（脚本 t0 或数字化轴解读问题），与 GEOS 无关；
> GEOS 应以 “validated analytical”（绿线）为基准比对。
