# DPDP Mandel 试错历史与文件清理记录

本目录清理后，只保留两份同源验证输入：

| 案例 | 输入文件 | 用途 |
|---|---|---|
| FIM 同源 | `DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml` | 全耦合法，同源 pressure 标定初始压力和位移荷载，`crossStorageOffDiagScale="1.0"`。这是当前主对比中的 FIM 输入。 |
| Seq 同源 | `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml` | 迭代耦合法，同源 pressure 标定初始压力和位移荷载，`crossStorageOffDiagScale="1.0"`。这是当前主对比中的 Seq 输入。 |

已删除的历史三案例曾使用的输出目录：

- FIM 0.911：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/fim0911`
- FIM 无修正：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/fim1000`
- 当前 Seq：`/tmp/dpdp_mandel_final_check_after_seq_strategy_20260715/seq`

历史三案例 XML 和对比图已不再长期保留；当前主图见下面同源 pressure 标定验证。

2026-07-16 新增同源 pressure 标定验证：

- FIM 同源输入：`DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml`
- Seq 同源输入：`DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml`
- FIM 输出：`/tmp/dpdp_same_source_unscaled_final_after_rebuild_20260716_111346/fim`
- Seq 输出：`/tmp/dpdp_same_source_unscaled_final_after_rebuild_20260716_111346/seq`
- 英文图：`analitical_result/GEOS_pressure_sameSource_FIM_Seq_analytical_20260716_1114_EN.png`
- 中文图：`analitical_result/GEOS_pressure_sameSource_FIM_Seq_analytical_20260716_1114_CN.png`
- 新旧位移荷载对比：`analitical_result/load_curve_old_vs_sameSourcePressure_20260716_0313.png`

本轮只比较当前解析脚本、FIM 同源、Seq 同源，不再把 0.911 拟合输入放入主图。图例中不强调
`crossStorageOffDiagScale`，中文图将 FIM/Seq 分别标为“全耦合法”和“迭代耦合法”；输入文件仍是
`crossStorageOffDiagScale="1.0"` 的同源对照。

2026-07-16 最终重跑时曾出现 Seq 首步不收敛。核查发现失败使用的是 2026-07-13 的旧
`build/bin/geosx`，而 Seq 修复提交在其后；重新编译后，同一输入可完整跑通。此次还修正了
同源 XML 生成物的格式问题：生成注释必须放在 XML declaration 之后，`loadFunction0000000`
片段改为 `.txt`，避免被全仓库 XML validation 当作完整 XML 校验。

## 关键结果

同源 FIM 完成完整验证，为 297 个时间步、0 次切步。同源 Seq 完成完整验证，为 422 个时间步、
0 次切步。

最新对比中的归一化压力峰值：

| 案例 | Matrix 峰值 `p/p0+` | Fracture 峰值 `p/p0+` |
|---|---:|---:|
| FIM 0.911 | 1.0484348 | 1.0376764 |
| FIM 无修正 | 1.0500023 | 1.0396572 |
| 当前 Seq | 1.0444486 | 1.0273431 |

此前 `final_pressure_compare.png` 和 `GEOS_pressure_compare_FIM0911_FIM1000_Seq_analytical_20260715_0928.png`
使用的是 `script/dpdp_mandel_analytical.py` 曲线。2026-07-15 后续检查发现，该脚本曲线与手动
Fig. 5c CSV 不一致，因此下面改用 workflow 指定的手动 CSV 作为压力主基准。

按手动 Fig. 5c CSV 统计的误差：

| 案例 | Matrix 平均/最大误差 | Fracture 平均/最大误差 |
|---|---:|---:|
| FIM 0.911 | 0.0542 / 0.3426 | 0.0171 / 0.0818 |
| FIM 无修正 | 0.0900 / 0.3056 | 0.0152 / 0.0728 |
| 当前 Seq | 0.0566 / 0.4484 | 0.0683 / 0.3401 |

采样点证据：

| `tau` | 手动 matrix | FIM 0.911 | FIM 1.0 | Seq | 手动 fracture | FIM 0.911 | FIM 1.0 | Seq |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 1.0397 | 1.0267 | 1.0150 | 1.0441 | 0.8218 | 0.8464 | 0.8408 | 0.9841 |
| 0.3 | 0.9583 | 0.9494 | 0.9080 | 1.0003 | 0.3169 | 0.3957 | 0.3867 | 0.6565 |
| 1.0 | 0.9094 | 0.8864 | 0.8238 | 0.8995 | 0.0155 | 0.0228 | 0.0215 | 0.1199 |
| 1000 | 0.7611 | 0.8763 | 0.8148 | 0.8707 | -0.0022 | 0.0000 | 0.0000 | 0.0000 |
| 3000 | 0.3262 | 0.6525 | 0.6078 | 0.7363 | -0.0027 | 0.0000 | 0.0000 | 0.0000 |

结论：历史三个 deck 都能完整跑通，但 FIM 和 Seq 都还没有严格对齐解析解。Seq fracture
排水偏慢仍然明确；matrix 晚期偏差目前受手动 CSV 与解析脚本不一致影响很大，需要先确认
最终压力基准，再继续判断 FIM/Seq 的物理误差来源。

2026-07-16 追查发现，历史 deck 的平台期偏低不是 FIM/Seq 程序必须使用 `0.911` 的证据，
而是初始压力和位移荷载混源造成的。历史 XML 使用 `p_m0/p_f0=0.455/0.488 Pc`，但旧位移
荷载初值 `|u_z0|=1.796142e-05 m` 对应 `2G|u_z0|/(P_c b)≈0.372`；按当前同源解析脚本与该
压力基准闭合，应为 `|u_z0|≈2.50037e-05 m`、`2G|u_z0|/(P_c b)≈0.518`。新旧荷载表时间点相同，
新表基本等于旧表整体乘以约 `1.392`，归一化形状差异约 `1e-4`。

同源 pressure 标定后的代表采样点：

| `tau` | 解析 matrix | FIM 同源 | Seq 同源 | 解析 fracture | FIM 同源 | Seq 同源 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 1.0597 | 1.0535 | 1.0549 | 0.9053 | 0.8766 | 0.8905 |
| 0.3 | 0.9668 | 0.9654 | 0.9667 | 0.4297 | 0.4225 | 0.4223 |
| 1.0 | 0.8882 | 0.8867 | 0.8942 | 0.0313 | 0.0265 | 0.0296 |
| 3.0 | 0.8827 | 0.8820 | 0.8896 | 0.0001 | 0.0002 | 0.0001 |
| 10.0 | 0.8840 | 0.8833 | 0.8908 | 0.0001 | 0.0000 | 0.0000 |
| 100.0 | 0.8905 | 0.8898 | 0.8977 | 0.0000 | 0.0000 | 0.0000 |

结论：FIM 同源平台期已与当前解析脚本一致；Seq 同源稳定跑通且平台期明显改善，
但 matrix 平台期相对 FIM/解析略高，晚期排水仍有差异。后续 Seq 问题应在同源荷载基准下继续查，
不应再用 `crossStorageOffDiagScale=0.911` 或其他非物理缩放补偿。

PDF 直接核查：用 PyMuPDF 渲染论文第 12 页后，Fig. 5c primary 曲线确实在横轴约
`10^3-10^4` 进入快速下降，手动 CSV 的晚期下降趋势与论文图一致。当前解析脚本和 GEOS
曲线更晚下降，因此下一步应优先修正解析脚本/位移加载来源，或改用独立应力加载验证。

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

`analitical_result/` 只长期保留手动解析 CSV、当前主对比图和用于说明历史混源错误的 load 曲线。
2026-07-15 清理掉该目录中旧阶段 PNG，包括旧 FIM/Seq 单案例图、digitized 对比图、BROKEN/历史
kappa 试验图和旧 Fig5c 诊断图。2026-07-16 继续删除旧 total-stress/scale 对比图、旧三案例 Fig5c
图、旧 sameSourcePressure 临时对比图，以及 `script/` 目录下早期 `GEOS_vs_analytical*.png` 截图。
当前保留图片为：

- `analitical_result/GEOS_pressure_sameSource_FIM_Seq_analytical_20260716_1114_CN.png`
- `analitical_result/GEOS_pressure_sameSource_FIM_Seq_analytical_20260716_1114_EN.png`
- `analitical_result/load_curve_old_vs_sameSourcePressure_20260716_0313.png`

## 已清理的历史输入表

`mandel_input_tables/` 保留基础初始位移表，以及同源 pressure 验证 XML 引用的
`same_source_pressure_calibrated/` 表：

- `xlin.geos`, `ylin.geos`, `zlin.geos`, `ux.geos`
- `xlin2.geos`, `ylin2.geos`, `zlin2.geos`, `uz.geos`

2026-07-15 删除未被当前保留 XML 或脚本引用的旧加载/诊断表：

2026-07-16 删除 `same_source_displacement_calibrated/`，该目录只是检查旧位移标定来源的中间生成物；
历史错误已记录在本文和 findings 中，最终同源 XML 只引用 `same_source_pressure_calibrated/`。

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
4. 当前保留的同源 FIM/Seq 都使用 `crossStorageOffDiagScale="1.0"`；历史 0.911 拟合输入已删除。
5. 2026-07-15 发现手动 Fig. 5c CSV 与解析脚本曲线不一致；继续调求解器前应先确认最终基准。
