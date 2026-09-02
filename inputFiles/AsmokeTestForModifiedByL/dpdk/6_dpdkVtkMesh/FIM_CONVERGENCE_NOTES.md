# Compositional Dual-Continuum Poromechanics FIM Convergence Notes

本文记录 `compositionalDualContinuumPoromechanics_vtk.xml` 的 FIM 收敛问题、代码修复、当前程序细节和输入文件配置方式。

## 1. 问题现象

目标算例：

```text
inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkVtkMesh/compositionalDualContinuumPoromechanics_vtk.xml
```

这是一个双网格双重介质模型：

- `mesh1/matrixRegion`: matrix continuum，同时承载力学自由度。
- `mesh2/fractureRegion`: fracture continuum，只有流动自由度。
- `matrixFlow` 和 `fractureFlow`: 两个 `CompositionalMultiphaseFVM`。
- `dualFlow`: matrix/fracture compositional cross-flow。
- `poroSolver`: `CompositionalMultiphaseDualContinuumPoromechanics`，使用 `FullyImplicit` 耦合。

原始问题分两类：

1. 早期严格 FIM 下 fracture volume residual 被污染，`fracture Rvol` 可停在 `O(1e-2)`。
2. 修掉明显污染后，仍存在 compositional volume constraint 尾部两点振荡，典型表现为：

```text
(fracture Rvol) ~= 3.15e-3 / 1.75e-6
(Rsolid)        ~= 7.77e-3
```

这导致 10 s 大步会到 `newtonMaxIter=100` 后切步。旧状态下完整 1000 s 算例虽然可以通过切步跑完，但有：

```text
Time step cuts                  = 2
Discarded nonlinear iterations  = 200
Successful nonlinear iterations = 498
```

## 2. 根因分析

### 2.1 Volume row 不应参与跨单元/跨物理耦合

Compositional FVM 每个 cell 的自由度通常为：

```text
0              : pressure
1..numComp     : component density variables
numComp        : local phase-volume closure equation row (Rvol)
```

这里的 volume closure row 是局部约束，不应像 mass row 一样参与邻居 flux 或外部力学交叉耦合。GEOS 通过 `globallyCoupledComponents` 标记哪些 row 需要全局耦合；`CompositionalMultiphaseFVM::setupDofs()` 会对 `numComp` 对应的 volume row 调用 `disableGlobalCouplingForEquation(...)`。

旧的 dual-continuum 稀疏模式构造在部分路径中仍按 `field.numComponents` 全部 row 建耦合，导致 fracture volume row 可能获得不应存在的跨介质或力学列。

### 2.2 部分装配路径对缺失矩阵列使用不安全写入

旧实现中存在类似：

```cpp
localMatrix.addToRowBinarySearchUnsorted(...)
```

如果目标 row 的 sparsity pattern 中没有对应 column，Release 下会出现写不到预期列甚至污染邻近数据的风险。对 FIM 这类强耦合问题，少量 row/column 污染就足以表现为 volume residual 尾部锁死。

### 2.3 fracture compositional volume closure 在 Newton 更新后没有保持当前态一致

`CompositionalMultiphaseBase::updateState()` 会根据当前 pressure/component density 更新 flash 变量和 phase volume fraction。FIM line search 和 Newton 每次更新 solution 后都会调用 `updateState()`。

问题是 fracture continuum 的 composition/phase volume closure 在 FIM 过程中会出现小偏离。如果只在初始化或收敛保存时闭合，而不在每次当前态更新后闭合，Newton residual 评估会在两个不一致状态之间来回跳：

- 一个状态接近 volume closure；
- 另一个状态因 mechanics/flow 更新后 fracture phase volume sum 偏离，重新激发 `Rvol`。

这就是后期 `Rvol` 两点振荡的直接原因。

## 3. 代码修复

### 3.1 DofManager: dual-continuum 稀疏模式排除 volume row

文件：

```text
src/coreComponents/linearAlgebra/DofManager.cpp
```

修改点：

- `setSparsityPatternDualContinuum(...)` 中 row loop 改为遍历 `globallyCoupledComponents`，不再遍历全部 `field.numComponents`。
- `countRowLengthsDualContinuum(...)` 同样只为 `globallyCoupledComponents` 统计 row length。
- `setSparsityPatternDualContinuumMechanics(...)` 和 `countRowLengthsDualContinuumMechanics(...)` 中，fracture flow row/col 也只使用 `fracPField.globallyCoupledComponents`。

效果：

- matrix/fracture cross-flow 稀疏模式不再包含 compositional volume row。
- mechanics-fracture flow block 不再把 displacement column 接到 fracture volume row。
- K_upf/K_pfu 只落到真实质量方程相关 row。

### 3.2 Cross-flow kernel: 只向已存在矩阵列写入

文件：

```text
src/coreComponents/physicsSolvers/multiphysics/dualContinuumCrossFlow/kernels/compositionalMultiPhase/FluxComputeKernel.hpp
```

新增 helper：

```cpp
void addToExistingRowEntries( localIndex localRow,
                              globalIndex const * colsToAdd,
                              real64 const * valsToAdd,
                              localIndex numCols ) const
```

该函数扫描 `m_localMatrix.getColumns(localRow)`，只对 row 中已经存在的 column 执行 `atomicAdd`。

替换位置：

- matrix side cross-flow mass rows。
- fracture side cross-flow mass rows。

效果：

- cross-flow Jacobian 不再依赖 `addToRowBinarySearchUnsorted` 对缺失列的行为。
- 若某 column 因 volume row 局部化或稀疏模式限制不存在，则安全跳过。

### 3.3 K_pfu: fracture mass row 对 displacement 的线性化补齐并安全写入

文件：

```text
src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp
```

关键函数：

```cpp
assembleFractureToMechanicsCouplingMultiphase(...)
```

物理含义：

fracture porosity 通过 matrix volumetric strain 更新：

```text
phi_f = phi_f,n + alpha_f * dEps_v + storage_terms
```

因此 fracture compositional mass equation 对 displacement 有导数：

```text
dM_c / dU = compDens_c * alpha_f * gradN * detJxW
```

修复内容：

- fracture row DOF 使用 fracture subregion 上的 flow dof key，而不是 matrix wrapper 上的辅助 dof。
- 只写 component mass rows `0..numComp-1`。
- volume-balance row `numComp` 不写 K_pfu。
- 写矩阵时扫描已有 row column，只对存在列 `atomicAdd`。
- 对 `useTotalMassEquation` 做 row transform：row 0 对应 total mass，其余 row 按 component shift。

效果：

- K_pfu 与 fracture accumulation 的 row 定义一致。
- 避免缺失 column 污染 fracture `Rvol`。

### 3.4 移除重复/不一致的 fracture volume closure Jacobian 装配

旧代码在 secondary flow assembly 后额外调用：

```cpp
assembleFractureVolumeClosureJacobianMultiphase(...)
```

该路径已不再调用。原因：

- compositional accumulation kernel 本身已经装配 volume balance residual/Jacobian。
- 额外手动补 fracture volume closure Jacobian 有重复或不一致风险。
- 正确策略是保持状态变量满足 closure，而不是给 volume row 添加跨物理补丁。

### 3.5 FIM updateState 后闭合 fracture compositional volume

文件：

```text
src/coreComponents/physicsSolvers/multiphysics/dualContinuumPoromechanics/DualContinuumPoromechanicsSolverBase.hpp
```

当前逻辑：

```cpp
virtual void updateState( DomainPartition & domain ) override
{
  Base::updateState( domain );

  if constexpr ( isMultiphaseFlow )
  {
    if( this->getNonlinearSolverParameters().couplingType() ==
        NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      enforceFractureCompositionalVolumeClosure( domain, false );
    }
  }
}
```

`enforceFractureCompositionalVolumeClosure(domain, false)` 做的事：

1. 遍历 fracture regions。
2. 计算每个 cell 的 `sum(phaseVolumeFraction)`。
3. 若偏离 1，则按

```text
globalCompDensity_c <- globalCompDensity_c / sum(phaseVolumeFraction)
```

缩放所有 component density。
4. 调用 fracture flow solver 的 `updateFluidState(fracSR)` 重新 flash。
5. 最多迭代 6 次，直到 closure error 小于 `1e-10`。
6. `saveConvergedState=false`，所以不会修改 `*_n` 或 converged history。

这一步只保持当前 Newton/line-search 状态一致，不改变时间步起点。

### 3.6 初始化时保存 fracture closure converged state

在 `initializePostInitialConditionsPostSubGroups()` 中，对于 multiphase + FullyImplicit，仍调用：

```cpp
enforceFractureCompositionalVolumeClosure( domain, true );
```

这用于保证初始 fracture composition/phase state 是闭合的，并保存为初始 converged state。

## 4. 当前 FIM 装配流程

`CompositionalMultiphaseDualContinuumPoromechanics` 的 FIM 装配主序列为：

1. 更新 fracture porosity fixed-stress 项。
2. 装配 mechanics block。
3. 装配 fracture-mechanics coupling `K_upf`。
4. 装配 mechanics -> fracture mass coupling `K_pfu`。
5. 装配 matrix face flux。
6. 装配 fracture flow system。
7. 装配 matrix/fracture cross-flow coupling `K_pmpf/K_pfpm`。
8. 装配 multiphase cross-storage 项。

Newton solution scaling：

```text
scale = Base scaling * PhysicsSolverBase oscillation scaling * fimNewtonRelaxation
```

FullyImplicit 下 `fimNewtonRelaxation=0.5` 用于抑制 dual-porosity pressure block 的 period-2 模式。

每次 Newton 或 line-search 应用 solution 后：

1. `Base::updateState(domain)` 更新 porosity/permeability/fluid state。
2. 对 fracture compositional variables 做当前态 volume closure。

## 5. 输入文件配置

目标 XML 中关键配置如下：

```xml
<CompositionalMultiphaseFVM name="matrixFlow"
  discretization="fluidTPFA"
  temperature="368.15"
  useMass="1"
  targetRegions="{ mesh1/matrixRegion }">
  <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="40" lineSearchAction="None"/>
  <LinearSolverParameters directParallel="0"/>
</CompositionalMultiphaseFVM>

<CompositionalMultiphaseFVM name="fractureFlow"
  discretization="fluidTPFA"
  temperature="368.15"
  useMass="1"
  targetRegions="{ mesh2/fractureRegion }">
  <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="40" lineSearchAction="None"/>
  <LinearSolverParameters directParallel="0"/>
</CompositionalMultiphaseFVM>

<CompositionalMultiPhaseDualContinuumFVM name="dualFlow"
  flowSolverName="matrixFlow"
  flowSolverName1="fractureFlow"
  matrixRegionList="{ matrixRegion }"
  fractureRegionList="{ fractureRegion }"
  targetRegions="{ mesh1/matrixRegion, mesh2/fractureRegion }">
  <DualContinuumCrossFlow
    fractureSpacingLx="10.0"
    fractureSpacingLy="10.0"
    fractureSpacingLz="10.0"
    fractureVolumeFraction="0.1"
    matrixRegionList="{ matrixRegion }"
    fractureRegionList="{ fractureRegion }"/>
</CompositionalMultiPhaseDualContinuumFVM>

<CompositionalMultiphaseDualContinuumPoromechanics name="poroSolver"
  dualcontinuumSolverName="dualFlow"
  solidSolverName="solidMech"
  targetRegions="{ mesh1/matrixRegion, mesh2/fractureRegion }"
  fractureVolumeFraction="0.1"
  fimNewtonRelaxation="0.5"
  enableFractureMechanicsCoupling="1"
  enableFimCrossStorage="1"
  stabilizationType="None"
  isThermal="0">
  <NonlinearSolverParameters
    couplingType="FullyImplicit"
    newtonMaxIter="100"
    newtonMinIter="0"
    newtonTol="1.0e-6"
    maxAllowedResidualNorm="1e13"
    lineSearchAction="Attempt"
    lineSearchResidualFactor="2.0"
    maxTimeStepCuts="10"
    timeStepCutFactor="0.5"/>
  <LinearSolverParameters solverType="direct" directParallel="1"/>
</CompositionalMultiphaseDualContinuumPoromechanics>
```

### 5.1 推荐配置

必须项：

- `couplingType="FullyImplicit"`: 使用 monolithic FIM。
- `newtonMinIter="0"`: 允许 `NewtonIter 0` 已满足 `1e-6` 时直接接受时间步。
- `newtonTol="1.0e-6"`: 当前已验证可收敛。
- `fimNewtonRelaxation="0.5"`: 保持 FIM Newton under-relaxation。
- `enableFractureMechanicsCoupling="1"`: 启用 `K_upf/K_pfu`。
- `enableFimCrossStorage="1"`: 启用 multiphase cross-storage。
- `lineSearchAction="Attempt"`: 当前推荐，不要改成 `Require`。

### 5.2 为什么不用 `lineSearchAction="Require"`

测试发现 `Require` 会在第一个时间步反复 cut：

```text
NewtonIter 0: R ~= 3.61e-1
NewtonIter 1: R ~= 4.53
Line search @ 0.031: R ~= 4.41e-1
```

虽然 line search 后 residual 已明显下降，但仍大于 NewtonIter 0 的 residual，因此 `Require` 判定失败并切步。这个模型的 FIM Newton 路径需要允许早期 residual 临时增大，再通过 relaxation 和后续迭代下降，所以应使用 `Attempt`。

### 5.3 `newtonMinIter=0` 的作用

若默认 `newtonMinIter=1`，很多小步在 `NewtonIter 0` 已满足 `R < 1e-6` 仍会被强制继续迭代。继续迭代会触发 `fimNewtonRelaxation=0.5` 的更新，反而可能扰动已满足的状态。因此 FIM 算例使用：

```xml
newtonMinIter="0"
```

## 6. 验证结果

构建命令：

```bash
make -C build geosx -j8
```

验证命令：

```bash
build/bin/geosx \
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkVtkMesh/compositionalDualContinuumPoromechanics_vtk.xml
```

修复后完整 1000 s 运行结果：

```text
Time steps                     = 37
Time step cuts                 = 0
Successful nonlinear iterations = 469
Discarded nonlinear iterations  = 0
```

最后一个 100 s 时间步收敛到：

```text
(matrix Rmass Rvol)   = (1.33e-10, 7.08e-16)
(fracture Rmass Rvol) = (3.30e-12, 8.14e-16)
(Rsolid)              = 8.54e-07
(R)                   = 8.54e-07
```

对比修复前：

```text
Time step cuts                 = 2
Discarded nonlinear iterations = 200
```

说明 fracture volume residual 尾部振荡已消除，`1e-6` 精度下 FIM 可以无切步完成。

## 7. 调试判断方法

推荐查看日志中的以下项：

```bash
rg -n "Time:|Cycle:|NewtonIter|Rmass Rvol|Rsolid|New dt|Time step cuts|Discarded nonlinear" run.log
```

正常状态：

- fracture `Rvol` 在收敛尾部应降到 `O(1e-15)`。
- 不应长期在 `O(1e-3)` 和 `O(1e-6)` 之间交替。
- `New dt = ...` 不应出现在该算例的 10 s/100 s 大步阶段。

异常状态：

- `fracture Rvol ~= 3e-3 / 1e-6` 两点振荡。
- `Rsolid` 卡在 `O(1e-3)` 到 `O(1e-2)`。
- `NewtonIter` 达到 99 后出现 `New dt = ...`。

## 8. 注意事项

- 不要重新给 fracture volume row 添加外部 mechanics/cross-flow column。volume row 是局部 closure row。
- 不要恢复 unsafe `addToRowBinarySearchUnsorted` 到可能缺列的 cross-flow 或 K_pfu 写入路径。
- `enforceFractureCompositionalVolumeClosure(domain, false)` 只应在当前态更新后使用，不应保存 converged state。
- 初始条件阶段可以使用 `saveConvergedState=true`，用于建立一致初始 fracture fluid state。
- `fimNewtonRelaxation=1.0` 已测试，收敛更差，不推荐。

