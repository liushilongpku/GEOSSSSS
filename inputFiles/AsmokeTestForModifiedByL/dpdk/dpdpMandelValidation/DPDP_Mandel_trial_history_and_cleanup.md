# DPDP Mandel 试错历史与文件清理记录

本目录清理后，只保留能复现当前 Mandel-Cryer 压力对比的三个输入文件：

| 案例 | 输入文件 | 用途 |
|---|---|---|
| FIM 0.911 | `DPDP_N2_dispdriven_fim_eff_direct_mesh10.xml` | 全隐式耦合，10x10，direct 求解器，有效介质力学参数，`crossStorageOffDiagScale="0.911"`。这是当前与解析压力曲线拟合最好的 FIM 输入。 |
| FIM 无修正 | `DPDP_N2_dispdriven_fim_eff_direct_mesh10_noCorrection.xml` | 与 FIM 0.911 输入相同，只把 `crossStorageOffDiagScale` 改为 `1.0`。这是不使用经验 cross-storage offdiag 缩放的 FIM 对照。 |
| 当前 Seq | `DPDP_N2_dispdriven_seq_eff.xml` | Sequential 耦合，10x10，有效介质力学参数，`crossStorageOffDiagScale="1.0"`，压力外迭代变化容差 `100 Pa`。这是当前修复后的 Seq 验证输入。 |

最新对比使用的输出目录：

- FIM 0.911：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/fim0911`
- FIM 无修正：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/fim1000`
- 当前 Seq：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/seq`

最新对比图：

- `/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/final_pressure_compare.png`
- `analitical_result/GEOS_pressure_compare_FIM0911_FIM1000_Seq_analytical_20260715_0928.png`

## 关键结果

两个 FIM 输入都完成了完整验证，均为 297 个时间步、0 次切步。当前 Seq 输入完成完整验证，
为 422 个时间步、0 次切步。

最新对比中的归一化压力峰值：

| 案例 | Matrix 峰值 `p/p0+` | Fracture 峰值 `p/p0+` |
|---|---:|---:|
| FIM 0.911 | 1.0484348 | 1.0376764 |
| FIM 无修正 | 1.0500023 | 1.0396572 |
| 当前 Seq | 1.0444486 | 1.0273431 |

结论：当前 Seq 已恢复 Mandel-Cryer 过冲，并且不需要 pressure relaxation 即可跑完整程。但 Seq
相对解析解仍偏慢，主要体现在 fracture 排水和 matrix 晚期衰减。FIM 0.911 仍是保留的 FIM
输入中与解析压力曲线拟合最好的一个。

## Seq 求解策略对照

旧的 Seq 保守策略主要是 storage / cross-storage split 修正前的兜底，不应在当前 deck 中默认恢复。
2026-07-15 做过以下对照：

| Seq 策略 | 输出目录 | 时间步 | nonlinear iteration | 切步 | 总时间 |
|---|---|---:|---:|---:|---:|
| 旧保守基线 | 历史日志 | 约 776 | 4974 | 0 | 约 35.5 s |
| 当前默认策略 | `/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/seq` | 422 | 2006 | 0 | 15.5 s |
| `100 Pa` + 早期细时间步 | `/tmp/dpdp_mandel_seq_balanced_dt_20260715` | 734 | 2996 | 0 | 31.0 s |

早期细时间步对照没有改善 fracture 排水 crossing：

| 指标 | 当前默认策略 | 早期细时间步 |
|---|---:|---:|
| `p_f/p0 = 0.5` crossing | 0.421253 | 0.421902 |
| `p_f/p0 = 0.1` crossing | 1.080309 | 1.109396 |

因此当前保留 `maxSequentialPressureChange=100 Pa` 和现有时间步分段；后续排水偏慢应作为模型/分裂精度问题处理。

## 已清理的历史图片

`analitical_result/` 只长期保留手动解析 CSV 和带时间戳的当前主对比图。2026-07-15 清理掉该目录
中旧阶段 PNG，包括旧 FIM/Seq 单案例图、digitized 对比图、BROKEN/历史 kappa 试验图和旧 Fig5c
诊断图。当前保留图片为：

- `analitical_result/GEOS_pressure_compare_FIM0911_FIM1000_Seq_analytical_20260715_0928.png`

## 已清理的历史输入表

`mandel_input_tables/` 只保留当前三个 XML 实际引用的初始位移表：

- `xlin.geos`, `ylin.geos`, `zlin.geos`, `ux.geos`
- `xlin2.geos`, `ylin2.geos`, `zlin2.geos`, `uz.geos`

2026-07-15 删除未被当前保留 XML 或脚本引用的旧加载/诊断表：

- `u_z_0.5.geos`
- `u_z_analitical.csv`
- `u_zmeter.geos`
- `u_zmeter.txt`
- `u_zmeter0.1x.csv`
- `u_ztime.geos`
- `u_ztime.txt`
- `uz_0.1.geos`

这些文件对应历史加载函数、缩放位移或人工诊断表；当前 XML 的实际加载函数已内嵌在各输入文件中，
初始位移只依赖上面保留的 8 个 `.geos` 表。

## 已删除 XML 的内容与删除原因

已删除的 XML 都是历史试错输入。下面保留它们的用途和结论，避免以后还需要打开旧 deck 才能知道
它们曾经验证过什么。

| 已删除输入 | 试错类别 | 不再保留的原因 |
|---|---|---|
| `DPDP_N2.xml` | 原始 N=2 Sequential 基础算例，内禀材料输入 | 已被位移驱动、有效介质参数的当前 Seq 输入替代。 |
| `DPDP_N2_constant_porosity.xml` | 常孔隙度诊断，用于隔离储量项贡献 | 仅是诊断输入，不属于最终解析解对比。 |
| `DPDP_N2_dispdriven.xml` | 早期位移驱动 Sequential 内禀材料输入 | 已被 corrected loading 和有效介质输入替代。 |
| `DPDP_N2_dispdriven_correctLF.xml` | 修正加载函数后的 Sequential 试错输入 | 中间阶段输入；最终 Seq 设置已进入保留输入。 |
| `DPDP_N2_dispdriven_correctLF_kappa.xml` | 修正加载函数并对裂缝渗透率做体积加权的试错输入 | 属于物理参数补偿试验，不作为最终验证输入。 |
| `DPDP_N2_dispdriven_fim.xml` | 早期 FIM 内禀材料、gmres 输入 | 已被 direct 求解器、有效介质参数的 FIM 输入替代。 |
| `DPDP_N2_dispdriven_fim_eff.xml` | 20x20 FIM 有效介质、gmres、scale 0.911 输入 | 已被更快的 10x10 direct FIM 输入替代。 |
| `DPDP_N2_dispdriven_fim_eff_direct.xml` | 20x20 FIM 有效介质、direct、scale 0.911 输入 | 已被 10x10 direct FIM 输入替代；当前比较不再需要常规保留 20x20。 |
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10_INTRINSIC.xml` | intrinsic 输入自动均匀化一致性检查 | 只是验证自动均匀化路径；最终比较保留 effective 输入。 |
| `DPDP_N2_dispdriven_fim_eff_direct_mesh10_stressload.xml` | FIM 应力加载 / RigidBoundary 变体 | 加载问题不同，不能直接并入当前位移驱动压力验证。 |
| `DPDP_N2_dispdriven_fim_eff_mesh40.xml` | 40x40 FIM 网格收敛检查 | 仅用于网格检查，成本高，不作为常规验证输入。 |
| `DPDP_N2_dispdriven_fim_intrinsic_kappa.xml` | FIM 内禀材料 + 裂缝渗透率体积加权 | 属于补偿试验，不作为最终验证输入。 |
| `DPDP_N2_identical_material_fim.xml` | FIM 双介质相同材料退化检查 | 求解器诊断，不是 Mandel 解析解对比输入。 |
| `DPDP_N2_identical_material_seq.xml` | Seq 双介质相同材料退化检查 | 求解器诊断，不是 Mandel 解析解对比输入。 |
| `DPDP_N2_only_drainage.xml` | 排水过程专项 Sequential 试验 | 诊断范围太窄；最终验证使用完整 Mandel-Cryer 过程。 |
| `DPDP_N2_stressload.xml` | 早期应力加载 Sequential 输入 | 边界/加载设置不同，不保留为位移驱动验证输入。 |
| `DPDP_N2_stressload_fixed.xml` | 修正后的应力加载 Sequential 变体 | 边界/加载设置不同，不保留为位移驱动验证输入。 |
| `DPDP_N2_symmetric_degradation.xml` | 双介质对称退化试验 | 求解器诊断，不属于最终解析解对比。 |
| `RigidBoundary_elastic_test.xml` | 纯弹性 RigidBoundary 边界测试 | 边界条件诊断输入，不是当前 Mandel 验证输入。 |
| `RigidBoundary_hetero_test.xml` | 非均质纯弹性 RigidBoundary 边界测试 | 边界条件诊断输入，不是当前 Mandel 验证输入。 |

## 已合并的非 XML 试错记录

| 原文件 | 试错类别 | 结论与处理 |
|---|---|---|
| `AUTO_SCHUR_CROSS_STORAGE_PROGRESS.md` | 自动 Schur cross-storage 标定尝试 | 2026-06-11 尝试用 FIM 装配矩阵估计 `K_pu K_uu^{-1} K_up`，以自动替代 XML 中的 `crossStorageOffDiagScale="0.911"`。Mandel-like 模式得到的 corrected offdiag 约 `-5.08671e-10`，几乎等于未缩放值，未能达到目标 `-4.635e-10`。说明该单标量/单 rank Schur 估计没有捕捉当前 Q1/Mandel 离散下的缺失抵消；代码已回退，当前不启用自动 Schur correction。该历史记录已合并到本文，原独立文件删除，避免被误认为当前方案。 |

## 这轮排查留下的判断

1. 原始 Seq 卡住不是因为总时长、网格数或线性求解器本身。受控扫描显示，根因是 Sequential
   split 中把 cross-storage offdiag 项隐式放入 flow solve 后导致分裂迭代发散。
2. 当前 Seq 在 fixed-stress Sequential 模式下滞后 cross-storage offdiag 项；这是分裂稳定化，不是
   改物理系数。只要外迭代收敛，固定点在压力外迭代容差内不变。
3. 当前保留 Seq 输入不需要 pressure relaxation。
4. 无修正 FIM 保留 `crossStorageOffDiagScale="1.0"`，用于检查不使用经验缩放时的结果。
5. FIM 0.911 作为最佳拟合输入保留，用于和无修正 FIM、当前 Seq、解析解放在同一张图里比较。
