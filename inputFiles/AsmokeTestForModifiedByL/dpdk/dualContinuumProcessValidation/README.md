<!-- Purpose: index dual-continuum validation cases and state their implementation status. -->

# Dual-continuum process validation

本目录保存双重介质交叉流基础机制、Thomas 1983 单块问题和 SPE6 重力排驱的输入、参考数据、脚本、
结果与报告。当前总体状态为：**case1-case7 已实现，case8 尚未实现。**

## Case status

| 编号 | 目录 | 内容 | 状态 | 报告 |
| --- | --- | --- | --- | --- |
| 1 | `1_pressure_exchange_combined` | 等体积/不等体积压差交换 | 已实现，通过 | [报告](1_pressure_exchange_combined/report_zh.md) |
| 2 | `2_gravity_exchange` | 单相重力驱替交换 | 已实现，通过 | [报告](2_gravity_exchange/report.md) |
| 3 | `3_capillary_exchange` | 毛管渗吸动力学 | 已实现，通过 | [报告](3_capillary_exchange/report.md) |
| 3.1 | `3.1_analytical` | 毛管渗吸闭式解析基准 | 已实现，通过 | [报告](3.1_analytical/report_zh.md) |
| 4 | `4_potgrad_derivative` | `PotGrad` 有限差分导数检查 | 已实现，通过（非运行时） | [报告](4_potgrad_derivative/report.md) |
| 5 | `5_Thomas_water_oil` | Thomas 水油渗吸细网格与单胞 | 已实现，复现成功 | [细网格](5_Thomas_water_oil/fine_grid/report_zh.md) / [单胞](5_Thomas_water_oil/single_block/report_zh.md) |
| 6 | `6_Thomas_oil_gas_depletion` | Thomas 气油枯竭细网格与单胞 | 已实现，复现成功 | [细网格](6_Thomas_oil_gas_depletion/fine_grid/report_zh.md) / [单胞](6_Thomas_oil_gas_depletion/single_block/RESULTS.md) |
| 7 | `7_SPE6_fine_grid_gravity_drainage` | SPE6 气油重力排驱细网格 | 已实现，验证成功 | [报告](7_SPE6_fine_grid_gravity_drainage/report_zh.md) |
| 8 | `8_SPE6_single_block_drainage` | SPE6 双重介质单胞 | 尚未实现 | [探索记录](8_SPE6_single_block_drainage/report_zh.md) |

case8 目录中的输入和试算用于诊断，当前约 `23.04%` 的 5 年体积采收率尚未达到 SPE6 约 `40%`
参考值，因此不能视为完成或通过。

## Evidence standard

每个完成的案例应说明验证目的、目标、控制方程、输入参数、参考解或文献基准、运行命令、实际结果、
误差和适用范围。不同案例的证据强度不同：

- case1-case3 使用 GEOS 运行历史和解析/守恒检查。
- case3.1 使用严格闭式解。
- case4 使用独立有限差分代数检查，不代表完整求解器运行验证。
- case5-case7 使用文献曲线、细网格/单胞交叉验证或独立 oracle。
- case8 尚未建立可接受的 canonical 实现和自动验收基线。

完整的案例说明、关键结果和限制见 [`summary.md`](summary.md)。

## Reproduction notes

case6 单胞提供统一复现入口：

```bash
python3 6_Thomas_oil_gas_depletion/single_block/scripts/reproduce.py
python3 6_Thomas_oil_gas_depletion/single_block/scripts/reproduce.py --run
```

运行产生的大型 VTK、HDF5、restart 和日志通常写入 `/tmp`，不作为验证输入提交。仓库保留能够重新生成
结果的 XML、表格、参考曲线、脚本、紧凑 CSV、图片和报告。
