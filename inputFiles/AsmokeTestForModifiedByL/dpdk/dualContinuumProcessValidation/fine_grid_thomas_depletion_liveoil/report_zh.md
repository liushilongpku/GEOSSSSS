<!-- 目的：汇总 Thomas 1983 10 ft 单重介质细网格验证。关键信息：最终采收率为 46.627%，运行原始数据不进入 Git。 -->

# Thomas 1983 10 ft 细网格验证总报告

## 1. 验证对象

本验证复现 Thomas（1983）论文中一个被含气裂缝包围的 10 ft 立方基质块。模型采用 `7 x 7 x 8`
单重介质细网格，不是一个网格内同时包含基质和裂缝的双重介质单块模型。

外部气体压力以 `0.75 psi/day` 下降。气体在重力作用下由上部侵入基质，较重的油向下排出；毛管压力、
气油密度差、三相相对渗透率和活油 PVT 共同控制采收率。验证基准为论文 Fig. 4 的 `3D model`
曲线及其约 46% 的平台值。

参考文献：`Thomas_1983_Fractured_Reservoir_Simulation.pdf`。

## 2. 最终结果

最终案例位于：

`continuous_pc_pvt_hydrostatic_thomas_stone/`

| 时间 | GEOS | Thomas Fig. 4 读图基准 | 差值 |
| ---: | ---: | ---: | ---: |
| 0.1 年 | 6.050% | 5.0% | +1.050 pp |
| 0.5 年 | 25.171% | 27.0% | -1.829 pp |
| 1.0 年 | 38.298% | 40.0% | -1.702 pp |
| 1.5 年 | 44.177% | 44.15% | +0.027 pp |
| 2.0 年 | 46.332% | 46.0% | +0.332 pp |
| 2.5 年 | **46.627%** | **46.0%** | **+0.627 pp** |

终值差异落在约正负 1 个百分点的论文读图不确定度内。早期 0.5--1.0 年仍有约 1.7--1.8 个百分点
差异，因此结论是：**终值和主要物理机制验证通过，瞬态曲线近似复现；不能宣称逐点精确匹配。**

## 3. 主要尝试和问题定位

| 阶段 | 主要变化 | 2.5 年附近结果 | 结论 |
| --- | --- | ---: | --- |
| 三段毛管压力缩放 | 按三个时间区间近似界面张力变化 | 57.09% | 排油过快，未通过 |
| 连续毛管压力分段 | 每 0.1 年更新毛管压力 | 约 55.9% | Pc 分段粗糙不是主因 |
| 加裂缝静水梯度 | 引入外部气体静水压力 | 约 43.3% | 静水梯度是重要因素 |
| 动态 PVT 静水梯度 | 气体密度随压力和 PVT 更新 | 43.699%--43.987% | 仍低于论文约 2 pp |
| 严格数值设置 | `newtonTol=1e-6`，`maxEventDt=5e4 s` | 43.797% | 排除时间步和容差为主因 |
| Thomas 式（25）修正 | 水油支路按正确变量解释 | **46.627%** | 最终验证结果 |

旧输入把水油油相相渗表预先反转后按油饱和度查询。在三相状态下，这不等价于 Thomas 式（25）中
按含水饱和度查询的 $k_{row}(S_w)$，从而低估了实际油相流动能力。最终输入在本算例近恒定
$S_w\approx0.20$ 的有效范围内采用 $k_{row}\approx1$，并保持 Thomas 气油支路不变。

独立审计得到：

- GEOS 体积加权 $k_{ro}=0.183703$；
- Thomas 式（25）体积加权 $k_{ro}=0.184178$；
- 二者比值为 `0.997420`。

## 4. 其他验证证据

### 采收率

![采收率与 Thomas Fig. 4 对比](continuous_pc_pvt_hydrostatic_thomas_stone/analysis/figures/recovery_comparison.png)

### 裂缝动态静水压力

![裂缝动态静水压力检查](continuous_pc_pvt_hydrostatic_thomas_stone/analysis/figures/fracture_hydrostatic_check.png)

2.5 年末端压力函数误差约为 $1.01\times10^{-3}$ Pa，25 个分段全链最大绝对压力误差约为
$0.362$ Pa。

### 纵向饱和度

![2.5 年纵向饱和度](continuous_pc_pvt_hydrostatic_thomas_stone/analysis/figures/final_vertical_saturation.png)

最终状态呈现上部气饱和度高、下部油饱和度高的空间分布，符合气油重力排驱的物理图像。

## 5. 可复现文件

- 基础输入：`thomas_10ft_7x7x8_liveoil_depletion.xml`
- 最终输入生成：`continuous_pc_pvt_hydrostatic_thomas_stone/generate_decks.py`
- 25 个分段输入：`continuous_pc_pvt_hydrostatic_thomas_stone/decks/right/`
- 顺序重启动运行：`continuous_pc_pvt_hydrostatic_thomas_stone/run_segments.py`
- 后处理：`continuous_pc_pvt_hydrostatic_thomas_stone/analysis/analyze_results.py`
- 相渗审计：`continuous_pc_pvt_hydrostatic_thomas_stone/analysis/audit_stone_relperm.py`
- 数值结果：`continuous_pc_pvt_hydrostatic_thomas_stone/analysis/*.csv`
- 详细报告：`continuous_pc_pvt_hydrostatic_thomas_stone/report_zh.md`

生成、运行和分析命令：

```bash
python3 inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/fine_grid_thomas_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone/generate_decks.py
python3 inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/fine_grid_thomas_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone/run_segments.py --variant right
python3 inputFiles/AsmokeTestForModifiedByL/dpdk/dualContinuumProcessValidation/fine_grid_thomas_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone/analysis/analyze_results.py
```

GEOS 运行生成的 `runs/`、HDF5、VTK、重启动文件和日志不进入 Git。Git 只保存输入、脚本、汇总
CSV、验证图片和报告。
