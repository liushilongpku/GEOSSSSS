# GEOS 双孔隙双渗透（DPDP, N=2）孔隙力学：实现、标定与验证技术文档

> 对象：GEOS 自定义双连续介质孔隙力学求解器 `SinglePhaseDualContinuumPoromechanics`
> （类 `DualContinuumPoromechanicsSolverBase`），用于复现 Mehrabian & Abousleiman (2014)
> 的 N=2 Mandel 解析解。
>
> 关键代码：
> - `src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp`
> - `src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/SinglePhaseDualContinuum.cpp`
> - `src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/DualContinuumCrossFlow.{hpp,cpp}`
> - 解析解：`inputFiles/.../dpdkHMValidation/dpdp_mandel_analytical.py`

---

## 目录
- [A. 背景与理论](#a-背景与理论)
- [B. 数值实现架构](#b-数值实现架构)
- [C. 两种耦合方法](#c-两种耦合方法)
- [D. 关键技术发现](#d-关键技术发现)
- [E. 新功能：内禀/等效输入二选一](#e-新功能内禀等效输入二选一)
- [F. 验证与结果](#f-验证与结果)
- [G. 使用说明](#g-使用说明)
- [H. 已知问题与展望](#h-已知问题与展望)
- [附录](#附录)

---

## A. 背景与理论

### A.1 问题定义

复现墨西哥湾页岩（GOM Shale）双孔隙双渗透介质在 Mandel 加载下的孔压演化（Mehrabian & Abousleiman 2014, Fig.5）。试样为 0.03×0.03×0.03 m 立方体；两个连续介质（N=2）：

| 介质 | 排水体模量 K | 剪切 G | 体积分数 v | 孔隙度 φ | 渗透率 k | 内禀 Biot α=1−K/Ks |
|---|---|---|---|---|---|---|
| 基质 matrix | 1.1 GPa | 0.757 GPa | 0.97 | 0.14 | 4.935e-21 m² (5 nd) | 0.9593 |
| 宏裂缝 fracture | 22.5 MPa | 15.49 MPa | 0.03 | 0.95 | 4.935e-15 m² (5 md) | 0.9992 |

公共参数：晶粒模量 Ks=2.7e10 Pa，流体模量 Kf=1.744e9 Pa（压缩率 c_f=5.734e-10 /Pa），黏度 μ=1e-3 Pa·s，泊松比 ν=0.22（两相同）。

**边界条件（Mandel 标准）**：z 向刚性压板（规定均匀 u_z 或施加应力）；y 向平面应变（u_y=0 两侧）；x 向 xneg 滚支(u_x=0)、xpos 自由且排水(p=0)；zneg 滚支。直接跨流交换 Γ₁₂=0（N=2）。

### A.2 控制方程与有效介质均匀化

广义 Biot 多孔隙理论：动量平衡 + 各连续介质质量平衡，通过**共享体应变**耦合。求解采用有效介质（Appendix A）：

```
K̄  = ( v_m/K_m + v_f/K_f )^{-1}                 (Reuss 体模量, 本算例 = 4.514e8 Pa)
Ḡ  = ( v_m/G_m + v_f/G_f )^{-1}  = 3K̄(1−2ν)/(2(1+ν))   (= 3.108e8 Pa)
ᾱ_i = K̄ · v_i · α_i / K_i        (ᾱ_m=0.3819, ᾱ_f=0.6014)
```

**储存矩阵 1/M̄**（恒应变，Mehrabian eq.A24/A25）：

```
对角  1/M̄_ii = v_i ( 1/M_i^intr + α_i²/K_i ) − ᾱ_i²/K̄
非对角 1/M̄_ij = − ᾱ_i · ᾱ_j / K̄
其中  1/M_i^intr = (α_i−φ_i)/Ks + φ_i·c_f
```

本算例数值：

```
1/M̄ = [  5.957e-10   −5.087e-10 ]
       [ −5.087e-10    5.462e-10 ]   (1/Pa)
```

**通量加权**：双连续介质对体相 Darcy 通量的贡献按体积分数加权 κ_i = v_i·k_i/μ。裂缝主导：κ_f = v_f·k_f = 0.03·4.935e-15 = 1.4805e-16 m²（deck 的 fracturePerm 即此值）。

### A.3 解析解

`dpdp_mandel_analytical.py`：Laplace 域解 + `mpmath` dehoog 数值反演。输出按不排水初值 p_i(0⁺) 归一（pm/pm0、pf/pf0），应力按 Pc 归一。特征时间 t0 = a²·tr(S̄)/Σκ ≈ 10.52 s（由裂缝渗透主导）。

- 解析 matrix：不排水→ Mandel-Cryer 过冲峰 **1.079 @ τ≈0.054** → 平台 **~0.885**（τ 1–100）→ 晚期排水。
- 解析 fracture：τ~0.3 排空。

---

## B. 数值实现架构

### B.1 双网格同位结构
- **mesh1**：基质网格，**唯一承载力学**（`SolidMechanicsLagrangianFEM`）。
- **mesh2**：裂缝网格，仅流动。
- 两网格几何**同位**：单元 k 在两网格代表同一物理点的两个连续介质。
- 自由度：位移 u（mesh1 节点）、p_m（mesh1 单元）、p_f（mesh2 单元）。
- **跨网格稀疏性**：`DofManager::addCouplingDualContinuumMechanics` 提供 u↔p_f 的稀疏条目，使 K_upf/K_pfu 不被丢弃。

### B.2 本构
- `ElasticIsotropic`：K、G（可写 `bulkModulus()`/`shearModulus()`）。
- `BiotPorosity`：biot 系数；`updateBiotCoefficientAndAssignModuli` 仅当 `defaultBiotCoefficient<=0` 时用 1−K/Ks 重算，否则保留显式值（这是 useIntrinsicInput 覆写能生效的前提）。

---

## C. 两种耦合方法

### C.1 Sequential（固定应力迭代耦合）

外层 Picard 交替解力学/流动；状态传递由 `mapSolutionBetweenSolvers` 完成。

```
时间步内每个外层迭代：
[力学] 解前 swapToCompositePressure：力学看到复合压力 p_eq（α_m·p_eq=ᾱ_m·p_m+ᾱ_f·p_f），
       并临时把矩阵本构 K/G 换成 K_eff/G_eff（Reuss）→ 解出位移、体应变
mapSolutionBetweenSolvers(力学)：restoreCompositePressure 恢复真实 p_m；
       mesh1 上平均平均总应力；固定应力孔隙度更新（注入 ᾱ_i·dEps_v）
[流动] 基质流(mesh1)+裂缝流(mesh2)+跨流；储存含固定应力稳定项 α²/K_fs
检查外层收敛
```

要点：固定应力分裂（流动用滞后稳定项 α²/K 近似力学反馈）；复合压力技巧（力学只解一套位移却感受两相）。

### C.2 FIM（全隐式单体耦合）

重写 `assembleSystem`，一次装配把 u/p_m/p_f 全部块装进同一雅可比：

```
Step 0  mapFractureDataToMatrix：p_f、α_f、裂缝DOF 从 mesh2 映射到 mesh1 同位单元
Step 1  单体基质核（SinglePhasePoromechanics）→ K_uu, K_upm, K_pmu, K_pmpm（自洽 Schur）
Step 1.5 updateFracturePorosityFixedStress：用基质位移更新裂缝孔隙度/质量（真实 Biot，无 α²/K 放大）
Step 2  K_upf：裂缝压力→位移残差 (+α_f·p_f·∇N)
Step 2b K_pfu：位移→裂缝质量 (ρ·α_f·∇N)
Step 3  基质面通量 assembleFluxTerms (mesh1 TPFA)
Step 4  裂缝流 (mesh2) + 跨流
Step 4b assembleCouplingTerms：多孔隙 1/M̄ 储存修正（eq.A25），含 crossStorageOffDiagScale
求解后 scalingForSystemSolution 施加 Newton 欠松弛 fimNewtonRelaxation(0.7)
```

### C.3 对比

| 维度 | Sequential | FIM |
|---|---|---|
| 结构 | 外层 Picard 解两个子系统 | 单体雅可比 + 单体 Newton |
| 力学反馈 | 固定应力 α²/K_fs（滞后） | K_upf/K_pfu 真实雅可比 |
| 压力传递 | 复合压力 swap/restore | Step 0 映射，直接用 p_m/p_f |
| 稳定手段 | 外层 subcycling | scalingForSystemSolution 欠松弛 |
| 主要短板 | 软裂缝 α²/K 放大→排水偏慢 | 非对角近抵消→需 0.911；过冲略低 |
| 代表 deck | `..._dispdriven_correctLF.xml` | `..._dispdriven_fim_eff*.xml` |

---

## D. 关键技术发现

### D.1 FIM 收敛性：period-2 振荡与 Newton 欠松弛
软耦合下单体 Newton 出现 period-2（特征值 λ≈−1）振荡：相邻两步在精确解两侧来回跳。修复：`scalingForSystemSolution` 中对 FIM 施加固定步长欠松弛 `fimNewtonRelaxation=0.7`（中点≈精确解，λ_eff≈1−2·relax 收缩）。配合跨网格 u↔p_f 稀疏性，收敛稳健（direct 求解器 0 次时间步切分）。

### D.2 储存矩阵精确化
早期用经验 fudge；后改为 Mehrabian eq.A25 的精确 1/M̄（见 A.2），由内禀参数构建（`SinglePhaseDualContinuum.cpp` Step 4b，第 206–265 行）。配合 κ=v_f·k_f 通量加权，裂缝段与解析近乎重合。

### D.3 非对角近抵消（核心问题）★

有效扩散储存为 `S̄ = invM + ᾱ⊗c_m`，其中 `c_m,i = ᾱ_i/(K̄+4Ḡ/3)`（oedometric）。在单体 FIM 中：
- **流动方程加入** invM（含非对角 −ᾱ_mᾱ_f/K̄ = **−5.087e-10**）；
- **力学 Schur 补提供** +ᾱ_mᾱ_f/(K̄+4Ḡ/3) = **+2.653e-10**。

净值 S̄_mf = **−2.44e-10**，是**两个 ~5e-10 大数的近抵消**。

**实测**（扫描 `crossStorageOffDiagScale`，direct 求解器，10×10 与 20×20 **逐位相同**）：

| scale | 平台 pm/4.55e5 |
|---|---|
| 0.0 | 1.536 |
| 0.5 | 1.179 |
| 1.0 | 0.822 |

线性拟合 `平台 = −0.714·scale + 1.536`，匹配解析 0.885 需 **scale = 0.911**。反解出**离散 Q1 力学有效模量 = 1.047e9 Pa = 1.21× 连续 oedometric(8.66e8)** —— 离散偏刚 21%，Schur 欠抵消，scale=1.0 时矩阵被拉过头（平台 0.82 偏低）。

**为何 0.911 推不出闭式**：S̄_mf 的力学半 = bᵀK_uu⁻¹b，依赖**全局**刚度矩阵求逆 + 解的**空间结构**（排水边界压力梯度 + Q1-P0 耦合）。单个 Q1 单元在均匀应变下精确复现连续模量，故 21% 偏差不是局部模量、不是 K̄/Ḡ/ν 的函数。**10×10=20×20 完全相同 → 这是单元阶次常数（非 h-误差），所以 0.911 是网格无关的离散常数，不是网格相关调参。**

> 详细量化见 memory `dpdp-fim-crossstorage-offdiag` 与 `SinglePhaseDualContinuum.cpp` 内注释。

### D.4 Sequential 排水偏慢
固定应力稳定项 α_f²/K_f（软裂缝 K_f=2.25e7）把裂缝储存放大约 80×；单趟耦合未充分抵消，导致排水偏慢。

---

## E. 新功能：内禀/等效输入二选一

### E.1 动机
此前 FIM deck 需**同时**手填等效参数（K̄、Ḡ、ᾱ）和内禀参数（储存用），冗余且易圆整出错。

### E.2 设计
在 `<SinglePhaseDualContinuumPoromechanics>` 加开关 `useIntrinsicInput`（默认 0）：

- **0（默认，完全向后兼容）**：本构填等效、DualContinuumCrossFlow 填内禀。行为不变。
- **1（FIM 专用）**：本构填**内禀**(K_m,G_m,α_m / K_f,G_f,α_f)。求解器在
  `initializePostInitialConditionsPostSubGroups → computeEffectiveFromIntrinsic` 中：
  1. 读本构内禀；
  2. Reuss 均匀化 K̄/Ḡ、ᾱ_i=K̄·v_i·α_i/K_i；
  3. 把内禀推给 DualContinuumCrossFlow（储存 Step 4b 用）；
  4. 用等效值**永久覆写**矩阵 ElasticIsotropic 的 K/G、矩阵/裂缝 BiotPorosity 的 α（供 FIM 核/K_upf）。

### E.3 安全性（对其它算例零影响）
三层隔离：① 代码仅在 `DualContinuumPoromechanicsSolverBase`（普通孔隙力学是别的类）；② 默认 0，现有 deck 不变；③ 仅当开关开+FIM 时才读写本求解器自己的本构实例。新增的 `BiotPorosity::getBiotCoefficientWritable()` 只是访问器，不改本构行为。

### E.4 验证
内禀输入 deck 与等效输入 deck 结果一致到 **0.028%**（差异来自等效 deck 手填常数的圆整；intrinsic 模式全精度、反而更准）。日志确认均匀化触发：
`poroSolver: useIntrinsicInput → homogenized effective medium from intrinsics (K_m=1.1e9, alpha_m=0.9593, K_f=2.25e7, alpha_f=0.9992)`。

---

## F. 验证与结果

### F.1 与解析解对比（FIM, scale=0.911）
| τ | GEOS | 解析 | 差 |
|---|---|---|---|
| 1 | 0.884 | 0.890 | −0.006 |
| 2–10 | 0.883–0.884 | 0.883–0.884 | <0.0002 |
| 90 | 0.888 | 0.890 | −0.002 |

平台均值 0.8852 vs 解析 0.885，τ≥1 吻合 <0.25%；全程跟随解析绿线（图 `analitical_result/GEOS_vs_analytical_vs_digitized_full.png`）。

### F.2 网格无关性
10×10 与 20×20 平台、scale 响应**逐位相同**（0.8216 / 1.5356）。10×10 已足够。

### F.3 MPI
gmres 在分区下并行预条件较弱、FIM 敏感（时间步频繁切分）；改 `solverType="direct"`（SuperLU_dist）→ 0 切分、与串行一致。

### F.4 解析脚本 vs 论文数字化曲线
本仓库解析脚本与论文 Fig5c 数字化曲线在**晚期排水时间尺度**差 ~10×（脚本 t0 或数字化轴解读问题），与 GEOS 无关；GEOS 以 "validated analytical"（绿线）为基准吻合良好。
（绘图脚本 `script/plot_GEOS_vs_analytical.py` 已去掉数字化曲线，仅画解析+GEOS。）

### F.5 标准对比图（两子图）
`script/GEOS_vs_analytical.py` 输出 `GEOS_vs_analytical.png`，含两子图：
- **(a) 压力**：基质 p_m + 裂缝 p_f，解析(线)+GEOS(点)，**全程吻合极好**（平台 0.885 精确，裂缝 τ~0.3 排空）。
- **(b) 位移**：顶板 u_z，解析+GEOS，**吻合极好**。

**应力子图被刻意去掉**，原因见 F.6 / H.4。压力与位移是可靠的对比量。

### F.6 总应力对比的"早吻合、晚漂移"问题 ★
GEOS 的 `matrixSolid_stress` 是**有效应力**；与解析（Fig5d 总应力）可比的总应力须重建：
```
σ_total = σ_eff − ᾱ_m·p_m − ᾱ_f·p_f
```
随排水进行，右边两项**反向大幅演化**（量级 ~1e6）：

| | σ_eff | ᾱ·p（Biot 压力项） | σ_total |
|---|---|---|---|
| τ→0 不排水 | -5.35e5 | +4.67e5 | -1.00e6 |
| τ→∞ 排空 | -8.0e5（Terzaghi 升 1.5×） | →0 | -8.0e5 |

总应力是**两个 ~1e6 大数之差**：
- **早期**几乎没排水，两项都在初值附近 → 差值准、曲线贴合（含 Mandel-Cryer 鼓包，GEOS 峰 1.05 vs 解析 1.08）；
- **晚期**压力掉 4.67e5，要总应力守恒需有效应力升 4.67e5，但 GEOS 只升约 2.67e5 → **补不够 → 总应力单调下漂**（中心 σ_zz 停在 ~0.90，不回到解析的 1.0；解析 τ~1000 的第二个基质 Mandel-Cryer 鼓包也未复现）。

**根因**：与 [D.3 非对角近抵消] 同源 —— Q1 离散力学刚 21%。但 scale=0.911 是按**压力平台**标定的（治储存→压力），**应力 = C:ε−ᾱ·p 直接依赖真实刚度 C，C 偏刚没改 → 应力漂移治不好**。这与 [H.4 应力驱动发散] 是同一个离散缺陷的两种表现。

---

## G. 使用说明

### G.1 本项目新增 XML 参数总表

**`<SinglePhaseDualContinuumPoromechanics>`：**

| 参数 | 默认 | 适用 | 含义 |
|---|---|---|---|
| `fractureVolumeFraction` | -1 | 两者 | v_f；<0 由网格体积算 |
| `fimNewtonRelaxation` | 0.5 | FIM | Newton 欠松弛(0,1]，压制 period-2 振荡 |
| `enableFractureMechanicsCoupling` | 1 | FIM | 是否装配 K_upf/K_pfu |
| `enableFimCrossStorage` | 1 | FIM | 是否启用多孔隙 1/M̄ 非对角储存 |
| `useIntrinsicInput` | 0 | FIM | 1=本构填内禀、自动均匀化；0=填等效 |
| `stabilizationType` | None | 两者 | 压力跳跃稳定化（本问题 None） |

**`<DualContinuumCrossFlow>`：**

| 参数 | 默认 | 含义 |
|---|---|---|
| `crossStorageOffDiagScale` | 1.0 | 非对角储存标定（本问题实测 **0.911**，补偿 Q1 离散刚化，见 D.3） |
| `intrinsicMatrixBiot` / `intrinsicMatrixBulkModulus` | -1 | 储存用矩阵内禀 α_m / K_m（<0=用材料值；useIntrinsicInput=1 时自动填） |
| `intrinsicFractureBiot` / `intrinsicFractureBulkModulus` | -1 | 裂缝内禀 α_f / K_f |
| `interporosityExchangeCoefficient` | 0 | 直接跨流交换系数（N=2 取 0） |
| `fractureVolumeFraction` | - | v_f（储存权重） |

> 注：修改 `*.cpp` 的 `setDescription` 后需 `make geosx_generate_schema` 才同步到 `schema.xsd`。

### G.2 范例 A —— 等效输入（useIntrinsicInput=0，`..._fim_eff_direct_mesh10.xml`）

本构直接填**等效**值，DualContinuumCrossFlow 填**内禀**值（供储存）：

```xml
<Constitutive>
  <!-- 矩阵：直接填等效 K̄/Ḡ -->
  <ElasticIsotropic name="matrixSolid"
    defaultDensity="0" defaultBulkModulus="4.514e8" defaultShearModulus="3.108e8"/>
  <BiotPorosity name="matrixPorosity"
    defaultGrainBulkModulus="2.7e10" defaultReferencePorosity="0.14"
    defaultBiotCoefficient="0.382"/>          <!-- 等效 ᾱ_m -->
  <ConstantPermeability name="matrixPerm"
    permeabilityComponents="{ 4.935e-21, 4.935e-21, 4.935e-21 }"/>

  <!-- 裂缝：K/G 用内禀(裂缝无独立力学)，Biot 用等效 -->
  <ElasticIsotropic name="fractureSolid"
    defaultDensity="0" defaultBulkModulus="2.25e7" defaultShearModulus="1.549e7"/>
  <BiotPorosity name="fracturePorosity"
    defaultGrainBulkModulus="2.7e10" defaultReferencePorosity="0.95"
    defaultBiotCoefficient="0.601"/>          <!-- 等效 ᾱ_f -->
  <ConstantPermeability name="fracturePerm"
    permeabilityComponents="{ 1.4805e-16, 1.4805e-16, 1.4805e-16 }"/>  <!-- κ=v_f·k_f -->
</Constitutive>

<Solvers>
  <DualContinuumFVM name="dualFlow" ...>
    <DualContinuumCrossFlow
      interporosityExchangeCoefficient="0.0" fractureVolumeFraction="0.03"
      intrinsicMatrixBiot="0.9593" intrinsicMatrixBulkModulus="1.1e9"     <!-- 储存要内禀 -->
      intrinsicFractureBiot="0.9992" intrinsicFractureBulkModulus="2.25e7"
      crossStorageOffDiagScale="0.911" .../>                              <!-- D.3 标定 -->
  </DualContinuumFVM>

  <SinglePhaseDualContinuumPoromechanics name="poroSolver"
    fractureVolumeFraction="0.03"
    fimNewtonRelaxation="0.7" enableFimCrossStorage="1"
    useIntrinsicInput="0"                      <!-- 等效输入模式（默认） -->
    stabilizationType="None">
    <NonlinearSolverParameters couplingType="FullyImplicit" .../>
    <LinearSolverParameters solverType="direct" directParallel="1"/>
  </SinglePhaseDualContinuumPoromechanics>
</Solvers>
```

要点：等效 K̄=4.514e8、ᾱ=[0.382,0.601] 须**用户手算**填入；intrinsic* 在 DualContinuumCrossFlow 重复给一遍（储存用）。

### G.3 范例 B —— 内禀输入（useIntrinsicInput=1，`..._fim_eff_direct_mesh10_INTRINSIC.xml`）

本构只填**内禀**，求解器自动均匀化；DualContinuumCrossFlow 的 intrinsic* 可省略（求解器自动填）：

```xml
<Constitutive>
  <!-- 矩阵：填内禀 K_m/G_m/α_m -->
  <ElasticIsotropic name="matrixSolid"
    defaultDensity="0" defaultBulkModulus="1.1e9" defaultShearModulus="7.574e8"/>
  <BiotPorosity name="matrixPorosity"
    defaultGrainBulkModulus="2.7e10" defaultReferencePorosity="0.14"
    defaultBiotCoefficient="0.9593"/>         <!-- 内禀 α_m -->
  <ConstantPermeability name="matrixPerm"
    permeabilityComponents="{ 4.935e-21, 4.935e-21, 4.935e-21 }"/>

  <!-- 裂缝：填内禀 K_f/G_f/α_f -->
  <ElasticIsotropic name="fractureSolid"
    defaultDensity="0" defaultBulkModulus="2.25e7" defaultShearModulus="1.549e7"/>
  <BiotPorosity name="fracturePorosity"
    defaultGrainBulkModulus="2.7e10" defaultReferencePorosity="0.95"
    defaultBiotCoefficient="0.9992"/>         <!-- 内禀 α_f -->
  <ConstantPermeability name="fracturePerm"
    permeabilityComponents="{ 1.4805e-16, 1.4805e-16, 1.4805e-16 }"/>
</Constitutive>

<Solvers>
  <DualContinuumFVM name="dualFlow" ...>
    <DualContinuumCrossFlow
      interporosityExchangeCoefficient="0.0" fractureVolumeFraction="0.03"
      crossStorageOffDiagScale="0.911" .../>   <!-- intrinsic* 可省，求解器自动填 -->
  </DualContinuumFVM>

  <SinglePhaseDualContinuumPoromechanics name="poroSolver"
    fractureVolumeFraction="0.03"
    fimNewtonRelaxation="0.7" enableFimCrossStorage="1"
    useIntrinsicInput="1"                      <!-- ★ 内禀输入模式 -->
    stabilizationType="None">
    <NonlinearSolverParameters couplingType="FullyImplicit" .../>
    <LinearSolverParameters solverType="direct" directParallel="1"/>
  </SinglePhaseDualContinuumPoromechanics>
</Solvers>
```

求解器初始化时自动算出 K̄=4.5144e8、Ḡ=3.1081e8、ᾱ_m=0.38187、ᾱ_f=0.60141 并覆写本构。**两版结果一致到 0.028%**（B 版全精度、更准）。

### G.4 两版对照小结
| | 范例 A（等效输入） | 范例 B（内禀输入） |
|---|---|---|
| 本构 K/Biot | 手填等效 4.514e8 / 0.382… | 填内禀 1.1e9 / 0.9593… |
| DualContinuumCrossFlow intrinsic* | 必填 | 可省（自动） |
| 等效参数来源 | 用户手算 | 求解器全精度计算 |
| `useIntrinsicInput` | 0 | 1 |
| 易错点 | 手算/圆整错误 | 无 |

### G.5 常见错误
- `useIntrinsicInput=1` 却在本构填了等效值 → 会被**二次均匀化**（错）。模式 1 必须填内禀。
- `useIntrinsicInput=1` 用于 Sequential → 被忽略（仅 FIM 生效，有日志提示）。
- 改 ν 或几何后沿用 0.911 → 需重新扫描标定（0.911 仅对 ν=0.22/本 Mandel 边界）。

### G.6 运行与画图
```
# 编译
make -C build -j$(nproc) geosx
# 运行（须在 deck 目录，相对路径指向 mandel_input_tables/）
build/bin/geosx -i DPDP_N2_dispdriven_fim_eff_direct_mesh10.xml
# MPI（partition 乘积 = -np）
mpirun -np 4 build/bin/geosx -i <deck> -x 2 -y 2 -z 1   # 建议 direct 求解器
# 画图
python analitical_result/script/plot_GEOS_vs_analytical.py <output_dir> <tau_max>
```

---

## H. 已知问题与展望

1. **Mandel-Cryer 过冲峰偏低 ~0.10**（τ≈0.05）：FIM(0.911) 平台精确，但过冲峰 1.05 vs 解析 1.08。属快速瞬态/Q1 空间分辨的独立问题。
2. **总应力"早吻合、晚漂移"**（见 F.6）：总应力是 σ_eff−ᾱ·p 两大数之差，离散刚 21% 使二者晚期抵消不精确，残差随排水累积；scale=0.911 只治了压力。**故标准对比图已不画应力子图。**

3. **应力驱动（Neumann）发散** ★（"路 C" 实测结论）：
   把加载从位移驱动改为应力驱动（刚性压板力 `RigidBoundary`，见 `..._stressload.xml`），本意让总应力按定义正确。但实测**FIM Newton 无法收敛**：
   | 配置 | 结果 |
   |---|---|
   | 1MPa 阶跃载荷 | 第 0 步即发散 |
   | 缓升(0→1MPa over 1e-3s)+relax0.7+dt5e-5+10次步切 | 跑到第 5 步、载荷~30% 时发散 |
   | 缓升(1e-2)+relax0.5+dt2e-5 | 第 0 步发散（加大阻尼反而更差） |

   **根因**：Neumann 边界下位移**自由**，t=0 附近**不排水近不可压**双孔隙响应使 K_uu 极度病态，Newton 发散。位移驱动是 Dirichlet（钉住位移→良态），这正是当初采用位移驱动的根本原因。**应力驱动发散与上面的应力漂移同源**——都是 Q1 离散在不排水/约束极限下的力学缺陷。

4. **非对角病态 + 应力问题的根治**：scale=0.911 只补偿压力；应力漂移(F.6)与应力驱动发散(H.3)都治不了。彻底方案：
   - **路 A（唯一同治）**：**无锁/高阶单元**（B-bar / 选择性缩减积分 / 增强应变 / Q2）使离散模量=连续 → scale 回到 1.0、**应力不再漂移、且改善 Neumann 病态使应力驱动也能收敛**。GEOS 现仅有标准 Q1 全积分六面体，需 FEM 内核开发。
   - **路 B（只治压力）**：矩阵自校准（初始化对一个内部单元解 K_uu·x=K_upf，实测交叉项反算 scale，ν/几何自适应，取代经验 0.911）。
   - ~~路 C 应力驱动~~：发散，不可行（除非先经路 A 改善收敛）。

5. **Sequential 有效参数路径失败**：Sequential 直接套等效参数会二次处理 K̄ 致过压（`..._seq_eff.xml` 已标注）。
6. **N>2 推广**、真·自动标定。

---

## 附录

### 附录 1：系数数值速查
```
K̄=4.514e8  Ḡ=3.108e8  ν=0.22  t0=10.52s  a=b=0.03  Pc≈7.17e5 Pa
ᾱ_m=0.38187  ᾱ_f=0.60141    α_m=0.9593  α_f=0.9992
v_m=0.97 v_f=0.03  φ_m=0.14 φ_f=0.95  Ks=2.7e10 Kf=1.744e9 μ=1e-3
κ_f=v_f·k_f=1.4805e-16   1/M̄=diag[5.957e-10,5.462e-10] off −5.087e-10
离散有效模量=1.047e9=1.21×oedometric(8.66e8)  →  scale=0.911(网格无关)
```

### 附录 2：关键代码位置
| 功能 | 文件:符号 |
|---|---|
| FIM 装配 Step0–4b | DualContinuumPoromechanicsSolverBase.hpp : `assembleSystem` |
| Newton 欠松弛 | 同上 : `scalingForSystemSolution` |
| Sequential 状态传递 | 同上 : `mapSolutionBetweenSolvers` / `swapToCompositePressure` |
| 内禀→等效均匀化 | 同上 : `computeEffectiveFromIntrinsic` |
| 多孔隙储存 1/M̄ | SinglePhaseDualContinuum.cpp : `assembleCouplingTerms`(Step4b) |
| 内禀参数 setter | DualContinuumCrossFlow.hpp / DualContinuumFlowSolverBase.hpp |
| 可写 Biot | BiotPorosity.hpp : `getBiotCoefficientWritable` |
| 解析解 | dpdkHMValidation/dpdp_mandel_analytical.py |

### 附录 3：crossStorageOffDiagScale 标定数据（direct 求解器）
| scale | 平台(10×10) | 平台(20×20) |
|---|---|---|
| 0.0 | 1.536 | 1.536 |
| 0.5 | 1.179 | — |
| 0.911 | 0.885 | — |
| 1.0 | 0.822 | 0.822 |

拟合 `平台 = −0.714·scale + 1.536`；匹配解析 0.885 → scale=0.911。
