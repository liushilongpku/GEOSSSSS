# DPDP Mandel 当前问题与结论状态

本文只记录技术结论：历史中发现了什么问题、是否已经解决、采用了什么办法、还剩什么问题。
运行步骤见 `DPDP_Mandel_validation_workflow.md`；输入文件清理记录见
`DPDP_Mandel_trial_history_and_cleanup.md`。

## 当前状态

截至 2026-07-16，当前只保留同源 pressure 标定验证输入：

- `DPDP_N2_dispdriven_fim_eff_direct_mesh10_sameSourcePressure.xml`
- `DPDP_N2_dispdriven_seq_eff_sameSourcePressure.xml`

最新完整对比显示：

- 同源 FIM、同源 Seq 均可跑完整个 `maxTime=308412 s` 验证。
- 当前 Seq 已恢复 Mandel-Cryer 过冲，不需要 pressure relaxation。
- 旧的保守 Seq 策略主要是旧 storage / cross-storage split 问题下的兜底；当前不应作为默认策略恢复。
- 历史三份混源输入已删除，仅保留文字记录；当前验证使用同源 FIM/Seq 输入。
- 压力基准存在不一致：手动 Fig. 5c CSV 和 `script/dpdp_mandel_analytical.py` 曲线在 fracture
  排水、matrix 晚期衰减上有明显差异。后续精度判断必须明确使用哪一个基准。
- 历史三份 XML 混用了不同来源的初始压力与位移荷载：`p_m0/p_f0=0.455/0.488 Pc`，
  但旧位移表对应 `2G|u_z0|/(P_c b)≈0.372`；按当前同源解析脚本应为
  `2G|u_z0|/(P_c b)≈0.518`。该混源错误导致 FIM/Seq 无缩放平台期偏低。
- 历史 FIM 0.911 是经验 `crossStorageOffDiagScale=0.911` 拟合输入，不能当作无修正物理输入；
  该 XML 已删除，不再作为当前验证输入。
- 历史 Seq 与解析解的差异更明显，主要表现为 fracture 排水和 matrix 晚期衰减偏慢。
- 使用同源 pressure 标定荷载后，FIM/Seq `crossStorageOffDiagScale=1.0` 均可完整跑通；
  FIM 平台期已与当前解析脚本一致，Seq 平台期明显改善；matrix 平台略高，晚期排水仍有差异。

## 问题与状态

| 问题 | 状态 | 证据/症状 | 处理办法 |
|---|---|---|---|
| 早期 Sequential 外迭代没有真正检查压力变化 | 已解决 | 历史日志中外迭代显示 `Max pressure change ... 0.000 Pa`，导致一步 staggered 就被接受 | 修正外迭代状态保存路径，使 Sequential 可以基于 flow pressure increment 做收敛判断 |
| `RigidBoundary` 早期投影式刚性压板导致收敛失败 | 已解决，但不再作为当前主验证路径 | 历史 stress-load deck 中 solid residual 平台化，投影会把 Newton iterate 推离平衡 | 将刚性约束改为系统内约束/惩罚形式；当前主验证采用位移驱动 deck，避免把压板边界问题混入 DPDP 验证 |
| 解析解实现和时间尺度曾经不一致 | 部分解决 | 历史脚本出现错误的 `t0` 和 storage 符号处理；早期注释也曾把 `tau` 与实际时间对应错 | 当前 `script/dpdp_mandel_analytical.py` 使用 Mehrabian-Abousleiman N=2 解析解，`t0≈10.5157 s`；但该脚本曲线仍与手动 Fig. 5c CSV 不一致，压力主误差暂按手动 CSV 统计 |
| 旧 displacement/load-function deck 与当前验证目标混杂 | 已解决 | 历史中存在 stress-load、correctLF、intrinsic/effective、kappa 补偿等多个 deck，容易混用 | 当前只保留两份同源验证输入；历史 deck 的用途和删除原因记录在清理总结中 |
| FIM 早期发散/不收敛 | 已解决 | 历史 FIM 在首个加载步或后续 Newton 迭代中出现发散/振荡 | 修正 FIM 相关组装、自由度映射和求解路径后，当前两个 FIM deck 均能完整跑完，0 次切步 |
| 双孔隙 storage matrix 处理缺失/不一致 | 已解决到当前验证可用状态 | 历史 GEOS matrix pressure 长期过高，不能复现解析解中的 matrix 平台和过冲 | 引入 dual-continuum cross-storage 处理；当前同源 FIM/Seq 均使用 `crossStorageOffDiagScale=1.0` |
| Seq full run 卡住/发散 | 已解决 | 受控扫描显示，当 Sequential split 中 cross-storage offdiag 项隐式进入 flow solve 时，外迭代会发散或极慢 | 当前 Seq 在 fixed-stress Sequential 模式下滞后 cross-storage offdiag 项；这是分裂稳定化，不是物理系数缩放 |
| Seq Mandel-Cryer 过冲缺失 | 已解决 | 历史 Seq 曲线曾缺少 matrix/fracture 早期压力过冲 | 修正 Seq 中孔隙度更新后的 fluid mass 同步、应变/压力映射以及 cross-storage split 后，过冲恢复 |
| Seq 求解策略过于保守 | 已解决 | 旧基线约 776 步、4974 次 nonlinear iteration、约 35.5 s；当前策略 422 步、2006 次 nonlinear iteration、约 15.5 s；恢复早期细时间步到 734 步也未改善 fracture 排水 crossing | 保留 `maxSequentialPressureChange=100 Pa` 与当前时间步分段；不恢复旧 `10 Pa` 容差和全程细时间步 |
| 压力基准不一致 | 未解决 | 手动 CSV 与解析脚本在 `tau=0.1` fracture 约为 `0.8218` vs `0.9059`，在 `tau=1000-3000` matrix 晚期衰减也明显不同 | 绘图脚本默认改为手动 Fig. 5c CSV，并可用 `--show-script-analytical` 叠加脚本曲线；继续改求解器前应先确认最终基准 |
| 历史位移荷载与初始压力混源 | 已定位，已生成同源对照 | 旧 XML 使用 `p_m0/p_f0=0.455/0.488 Pc`，但 `loadFunction0000000` 初值 `1.796142e-05 m` 对应的初始压力约 `0.3265/0.3506 Pc`；旧荷载整体比同源 pressure 标定荷载小约 `1/1.392` | 新增 `generate_same_source_p0_uz.py`，生成同源 `p0/uz/loadFunction`；新增同源 FIM/Seq 无缩放输入并完成验证 |
| 位移驱动加载函数依赖可疑解析脚本 | 部分解决 | 已确认历史 deck 的压力和位移荷载混源；同源 pressure 标定后 FIM 平台期与当前解析脚本一致。但当前解析脚本与论文手动 Fig. 5c CSV 的晚期差异仍未解决 | 短期使用同源 pressure 标定 deck 做程序一致性验证；论文级验证仍需修正/确认解析脚本绝对标定和 Fig. 5c 基准 |
| FIM 与解析解仍不完全一致 | 未完全解决 | 按手动 Fig. 5c CSV 统计，FIM 0.911 matrix/fracture 平均误差约 `0.0542/0.0171`，最大误差约 `0.3426/0.0818`；FIM 1.0 为 `0.0900/0.0152`，最大误差约 `0.3056/0.0728`。matrix 最大误差主要来自手动 CSV 的晚期下降 | 先核查解析基准；基准确认后再继续查 cross-storage/offdiag 理论系数和离散化 |
| Seq 与解析解仍不完全一致 | 未完全解决 | 按手动 Fig. 5c CSV 统计，Seq matrix/fracture 平均误差约 `0.0566/0.0683`，最大误差约 `0.4484/0.3401`；fracture 在 `tau=0.1-1` 排水明显偏慢 | 保留为后续精度问题；不得用 pressure relaxation 或无物理意义的参数缩放掩盖 |

## 当前保留结论

1. 当前 Seq 的稳定性问题不是网格数、总时长或线性求解器导致，根因是 Sequential split 下 cross-storage offdiag 的隐式处理方式。
2. 当前 Seq 使用的 offdiag 滞后是分裂稳定化；外迭代收敛后，固定点在压力外迭代容差内保持一致。
3. 当前验证不需要 pressure relaxation。
4. 历史 FIM 0.911 是经验拟合输入，已删除；当前保留输入均使用 `crossStorageOffDiagScale=1.0`。
5. 历史平台期偏低的直接原因是初始压力与位移荷载混源；`crossStorageOffDiagScale=0.911`
   只是补偿该混源错误的经验拟合。
6. 当前最先需要处理的是解析脚本与论文图的基准差异，而不是继续调 FIM/Seq 参数。
7. 后续如果新增案例，应先用通用绘图脚本与明确基准对比，再把新结论更新到本文，而不是写入 workflow。
