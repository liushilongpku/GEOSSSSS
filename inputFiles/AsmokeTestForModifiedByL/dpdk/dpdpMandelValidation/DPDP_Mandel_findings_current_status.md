# DPDP Mandel 当前问题与结论状态

本文只记录技术结论：历史中发现了什么问题、是否已经解决、采用了什么办法、还剩什么问题。
运行步骤见 `DPDP_Mandel_validation_workflow.md`；输入文件清理记录见
`DPDP_Mandel_trial_history_and_cleanup.md`。

## 当前状态

截至 2026-07-15，当前保留的三个验证输入为：

- `DPDP_N2_dispdriven_fim_eff_direct_mesh10.xml`
- `DPDP_N2_dispdriven_fim_eff_direct_mesh10_noCorrection.xml`
- `DPDP_N2_dispdriven_seq_eff.xml`

最新完整对比显示：

- FIM 0.911、FIM 1.0、Seq 均可跑完整个 `maxTime=308412 s` 验证。
- 当前 Seq 已恢复 Mandel-Cryer 过冲，不需要 pressure relaxation。
- 旧的保守 Seq 策略主要是旧 storage / cross-storage split 问题下的兜底；当前不应作为默认策略恢复。
- 当前 Seq 与解析解仍有差异，主要表现为 fracture 排水和 matrix 晚期衰减偏慢。
- FIM 0.911 是当前与解析压力曲线拟合最好的 FIM 输入；FIM 1.0 是无经验 cross-storage offdiag 缩放的对照输入。

## 问题与状态

| 问题 | 状态 | 证据/症状 | 处理办法 |
|---|---|---|---|
| 早期 Sequential 外迭代没有真正检查压力变化 | 已解决 | 历史日志中外迭代显示 `Max pressure change ... 0.000 Pa`，导致一步 staggered 就被接受 | 修正外迭代状态保存路径，使 Sequential 可以基于 flow pressure increment 做收敛判断 |
| `RigidBoundary` 早期投影式刚性压板导致收敛失败 | 已解决，但不再作为当前主验证路径 | 历史 stress-load deck 中 solid residual 平台化，投影会把 Newton iterate 推离平衡 | 将刚性约束改为系统内约束/惩罚形式；当前主验证采用位移驱动 deck，避免把压板边界问题混入 DPDP 验证 |
| 解析解实现和时间尺度曾经不一致 | 已解决 | 历史脚本出现错误的 `t0` 和 storage 符号处理；早期注释也曾把 `tau` 与实际时间对应错 | 当前 `script/dpdp_mandel_analytical.py` 使用 Mehrabian-Abousleiman N=2 解析解，`t0≈10.5157 s`；解析资源路径已改为相对当前验证目录 |
| 旧 displacement/load-function deck 与当前验证目标混杂 | 已解决 | 历史中存在 stress-load、correctLF、intrinsic/effective、kappa 补偿等多个 deck，容易混用 | 当前只保留三个验证输入；历史 deck 的用途和删除原因记录在清理总结中 |
| FIM 早期发散/不收敛 | 已解决 | 历史 FIM 在首个加载步或后续 Newton 迭代中出现发散/振荡 | 修正 FIM 相关组装、自由度映射和求解路径后，当前两个 FIM deck 均能完整跑完，0 次切步 |
| 双孔隙 storage matrix 处理缺失/不一致 | 已解决到当前验证可用状态 | 历史 GEOS matrix pressure 长期过高，不能复现解析解中的 matrix 平台和过冲 | 引入 dual-continuum cross-storage 处理；FIM 保留 0.911 拟合输入和 1.0 无修正对照输入 |
| Seq full run 卡住/发散 | 已解决 | 受控扫描显示，当 Sequential split 中 cross-storage offdiag 项隐式进入 flow solve 时，外迭代会发散或极慢 | 当前 Seq 在 fixed-stress Sequential 模式下滞后 cross-storage offdiag 项；这是分裂稳定化，不是物理系数缩放 |
| Seq Mandel-Cryer 过冲缺失 | 已解决 | 历史 Seq 曲线曾缺少 matrix/fracture 早期压力过冲 | 修正 Seq 中孔隙度更新后的 fluid mass 同步、应变/压力映射以及 cross-storage split 后，过冲恢复 |
| Seq 求解策略过于保守 | 已解决 | 旧基线约 776 步、4974 次 nonlinear iteration、约 35.5 s；当前策略 422 步、2006 次 nonlinear iteration、约 15.5 s；恢复早期细时间步到 734 步也未改善 fracture 排水 crossing | 保留 `maxSequentialPressureChange=100 Pa` 与当前时间步分段；不恢复旧 `10 Pa` 容差和全程细时间步 |
| Seq 与解析解仍不完全一致 | 未完全解决 | 当前 Seq matrix 峰值约 `1.04445`，fracture 峰值约 `1.02734`；排水和晚期衰减相对解析解偏慢 | 保留为后续精度问题；不得用 pressure relaxation 或无物理意义的参数缩放掩盖 |

## 当前保留结论

1. 当前 Seq 的稳定性问题不是网格数、总时长或线性求解器导致，根因是 Sequential split 下 cross-storage offdiag 的隐式处理方式。
2. 当前 Seq 使用的 offdiag 滞后是分裂稳定化；外迭代收敛后，固定点在压力外迭代容差内保持一致。
3. 当前验证不需要 pressure relaxation。
4. FIM 0.911 是经验拟合输入；FIM 1.0 是无经验缩放对照。两者都需要保留，不能互相替代。
5. 后续如果新增案例，应先用通用绘图脚本与解析解对比，再把新结论更新到本文，而不是写入 workflow。
