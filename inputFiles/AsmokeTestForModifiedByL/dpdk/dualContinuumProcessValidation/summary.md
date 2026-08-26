<!-- Purpose: summarize the purpose, target, justification, evidence, and status of every dual-continuum validation case. 关键信息：除 6、8 未验证上，其余案例均验证或复现成功。 -->

# 双重介质过程验证总表

本目录对双重介质基质—裂缝交叉流相关机制进行拆开验证：压差驱动、重力驱替、毛管/渗吸驱替、
`PotGrad` 导数契约，以及两个细网格文献复现（Thomas 活油枯竭、SPE6 重力排驱）。它们不是
一个混合算例，而是按机制拆开的证据链。

## 状态总览

| 编号 | 案例 | 过程 | 验证级别 | 当前状态 |
| --- | --- | --- | --- | --- |
| 1 | pressure exchange | 压差驱动交换（等体积 + 不等体积） | GEOS 运行 + 两储集体解析解 | **通过** |
| 2 | gravity exchange | 重力驱替交换 | GEOS 运行 + 重力驱替势差平衡 | **通过** |
| 3 | capillary exchange | 毛管/渗吸交换（动力学版） | GEOS 组装诊断 + 毛管状态 + 动态响应与守恒判据 | **通过** |
| 3.1 | analytical | 毛管渗吸解析基准版 | 严格两储集体闭式解 + 数值对照 | **通过** |
| 4 | potgrad | `PotGrad` 毛管导数契约 | 独立有限差分代数检查 | **通过，但非运行时验证** |
| 5 | Thomas fine grid | Thomas 1983 活油细网格枯竭 | 论文 Fig. 4 采收率曲线 | **复现成功** |
| 7 | SPE6 fine grid | SPE6 单块气油重力排驱 | 论文 Fig. 1 采收率曲线（7×7×8 / 10×10×10） | **验证成功** |
| 6 | Thomas single block | Thomas 单块驱替 | 目标基准尚未建立 | **未验证** |
| 8 | SPE6 single block | SPE6 单块（单重介质单网格）驱替 | 论文约 40% 采收率 | **未验证上** |

说明：编号 `6` 与 `8` 目前未通过验证。其余案例（1、2、3、3.1、4、5、7）均已通过对应的
验证或复现判据。各案例的子报告、输入、脚本、图片与运行归档保留在各自目录中。

## 1：压差驱动交换（pressure exchange）

### 验证目的

隔离最基本的基质—裂缝压力交换，确认交叉流项的方向、两侧等量反向装配、储集量守恒和后向
Euler 时间推进。

### 验证目标与影像

两个封闭、共位、单相单元从初始压差出发，应满足：高压侧压力下降、低压侧压力上升、压差单调
衰减；加权平均压力不变；首步衰减因子等于解析值。

基质初始高压、裂缝低压，流体经交叉通道从基质流向裂缝，两端最终趋近同一平衡压力
（本算例约 `1.5 MPa`）。方程退化为

```text
Cm dpm/dt = -A (pm - pf)
Cf dpf/dt = +A (pm - pf)
```

因此压差解析解与守恒量可直接写出。

### 参考解和结果

输入为 `1_pressure_exchange_combined/` 下的等体积/不等体积两个算例；详细推导、图片与数值见
[`1_pressure_exchange_combined/report_zh.md`](1_pressure_exchange_combined/report_zh.md)。

```text
GEOS exit code               0
numerical q                  0.8928571053
reference q                  0.8928571429
relative q error             4.21e-8
relative mean-pressure drift 8.81e-8
resolved-curve max error     4.39e-7
```

在可分辨区间内误差满足判据，因此 1 通过；它不验证毛管、重力或 compositional `PotGrad`。

## 2：重力驱替交换（gravity exchange）

### 验证目的与目标

验证 Kazemi 型重力驱替势差（GDP）的符号、大小及加入矩阵侧势差的方向。初始两侧等压，但
`rho_m=1100`、`rho_f=800 kg/m3`、`g_z=-9.81 m/s2`、`Lz=100 m`。目标平衡压差

```text
p_m - p_f = -GDP,  GDP = 147150 Pa
```

### 影像与依据

两个连续体各一个共位单元，普通单元中心重力水头相互抵消，初始压差为零，剩余驱动力只能是
GDP；稳态势差为零直接得到上式。

### 结果

```text
reference GDP                 147150.0000 Pa
final p_m - p_f              -147150.0803 Pa
relative error                5.46e-7
```

判据相对误差 `<1e-5`、绝对误差 `<1 Pa` 均满足，因此 2 通过（单相验证）。详见
[`2_gravity_exchange/report.md`](2_gravity_exchange/report.md)。

## 3：毛管/渗吸交换动力学版（capillary exchange）

### 验证目的与依据

验证饱和度相关毛管压力、零重力下的渗吸过程。采用非零驱动 + 相分数变化 + 组分守恒的过程判据，
证明模块未破坏物理；由于其物性（CO2-brine flash、Brooks-Corey）无法化为闭式解，故采用动态
历史与守恒判据而非独立解析参考。

### 结果

GEOS 正常完成（`10449` 个接受步、零切步），动态历史与守恒判据通过。详见
[`3_capillary_exchange/report.md`](3_capillary_exchange/report.md)。

## 3.1：毛管渗吸解析基准版（analytical）

### 验证目的与依据

在动力学版之上，另建一个可严格闭式求解的两储集体毛管渗吸问题：水不可压缩、气体密度与黏度
恒定、相对渗透率近乎恒定、基质毛管压力线性、裂缝毛管压力为零。基质水饱和度偏离平衡值的量
严格满足一阶线性 ODE，解析解为单指数。

### 结果

用平衡点 `S_eq=0.6`、基质初始 Pc `5.0e5 Pa` 与线性下降构造的闭式解与数值响应一致，因此
3.1 通过。相渗/毛管曲线与推导见
[`3.1_analytical/report_zh.md`](3.1_analytical/report_zh.md) 及对应图片。

## 4：`PotGrad` 毛管导数契约（potgrad）

### 验证目的与结果

独立有限差分代数检查压力、饱和度、重力与毛管项的组合 Jacobian 契约。该检查通过但**不是运行时
验证**：它不能替代 1/2/3 的运行时字段映射证明。详见
[`4_potgrad_derivative/report.md`](4_potgrad_derivative/report.md)。

## 5：Thomas 1983 活油细网格枯竭（Thomas fine grid）

### 验证对象与结果

复现 Thomas（1983）中被含气裂缝包围的 `10 ft` 基质块，采用 `7 x 7 x 8` 单重介质细网格。
外部气体压力以 `0.75 psi/day` 下降，气体由上部侵入、较重油向下排出。验证基准为论文 Fig. 4
`3D model` 曲线约 46% 平台值。

最终案例位于 `5_Thomas_fine_grid_depletion_liveoil/continuous_pc_pvt_hydrostatic_thomas_stone/`：

| 时间 | GEOS | Thomas Fig. 4 | 差值 |
| ---: | ---: | ---: | ---: |
| 1.0 年 | 38.298% | 40.0% | -1.702 pp |
| 1.5 年 | 44.177% | 44.15% | +0.027 pp |
| 2.5 年 | **46.627%** | **46.0%** | **+0.627 pp** |

终值落在论文读图不确定度内，主要机制与终值通过，早期 0.5--1.0 年约有 1.7--1.8 个百分点
差异，因此结论为“终值和主要物理机制验证通过，瞬态曲线近似复现”。详见
[`5_Thomas_fine_grid_depletion_liveoil/report_zh.md`](5_Thomas_fine_grid_depletion_liveoil/report_zh.md)。

## 7：SPE6 单块气油重力排驱细网格（SPE6 fine grid）

### 验证对象与结果

在 SPE6 单块物理约束下（4500 psig、5 年、活油 PVT、零裂缝毛管压力、固定裂缝薄壳）用细网格
复现气油重力排驱。两种分辨率均为 VTK-only、不设 `maxEventDt`：

| 子目录 | 总网格 | 内部基质 | 5 年体积采收率 | 相对 SPE6 40% |
| --- | ---: | ---: | ---: | ---: |
| `7x7x8` | 392 | 150 | 38.960% | -1.040 pp |
| `10x10x10` | 1000 | 512 | 39.751% | -0.249 pp |

两者都接近 SPE6 约 40% 的工程目标，且保持“上部富气、下部富油、水相基本不动”的重力驱替
特征。因无时间步上限，数值含可见时间离散误差（10×10×10 若限制 100000 s 则约 39.924%），
但结论仍为验证成功。详见
[`7_SPE6_fine_grid_gravity_drainage/report_zh.md`](7_SPE6_fine_grid_gravity_drainage/report_zh.md)。

## 6：Thomas 单块驱替（未验证）

对应目录 `6_Thomas_single_block_drainage`。当前目录为空，尚未建立可对比的论文基准，也无
GEOS 运行结果，因此标记为**未验证**。需要后续补充 Thomas 单块解析/文献基准或运行算例后再更新。

## 8：SPE6 单块（单重介质单网格）驱替（未验证）

对应目录 `8_SPE6_single_block_drainage`（承接原 `P0_thomas_single_block`）。该单网格、
单重介质算例尝试复现 SPE6 单块约 40% 采收率，但反复试错后仍未与论文结果对上，因此标记为
**未验证上**。原 `P0_thomas_single_block` 的全部输入、脚本、图片与试错记录已迁入此目录；
后续若更换物性/相渗/毛管处理需重跑并重新登记状态。

## 综合说明

1. 编号 1/2/3/3.1 的“通过”只对各自隔离的单一机制负责。
2. 编号 4 的通过是代数契约检查，不能替代运行时字段映射证明。
3. 编号 5/7 是文献细网格复现，结论为“终值与主要机制通过、瞬态近似复现”，不是逐点精确匹配。
4. 编号 6/8 未验证，需在补充基准或修正物理设置后重跑再更新本表。
5. 各案例的运行输出归档在其目录内，不依赖公共输出目录；`/tmp` 运行数据不纳入 Git。
