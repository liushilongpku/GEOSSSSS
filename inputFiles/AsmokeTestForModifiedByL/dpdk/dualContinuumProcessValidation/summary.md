<!-- Purpose: summarize the purpose, target, justification, evidence, and status of every dual-continuum validation case. -->

# 双重介质过程验证总表

本目录针对双重介质基质-裂缝交叉流中三个过程分别验证：压差驱动、重力驱替和毛管/渗吸驱替；同时加入一个局部 `PotGrad` 导数契约检查，以及一个流固耦合 Mandel 基准。它们不是一个混合算例，而是按机制拆开的证据链。

## 状态总览

| 案例 | 过程 | 验证级别 | 当前状态 |
| --- | --- | --- | --- |
| A0 | 压差驱动交换（等体积 + 不等体积） | GEOS 运行 + 两储集体解析解；考虑压力差浮点分辨率 | **通过** |
| G0 | 重力驱替交换 | GEOS 运行 + GDP 解析平衡 | **通过** |
| I0 | 毛管/渗吸交换 | GEOS 组装诊断 + 毛管状态输出 + 动态响应和守恒判据 | **通过** |
| J0 | `PotGrad` 毛管导数 | 独立有限差分代数检查 | **通过，但非运行时验证** |
| P1 | 双孔隙双渗 Mandel 流固耦合 | 文献解析解 + 现有 FIM/Seq 工作流 | **有条件验证** |

## A0：压差驱动

### 验证目的

隔离最基本的基质-裂缝压力交换，确认交叉流项的方向、矩阵/裂缝两侧的等量反向装配、储集量守恒和后向 Euler 时间推进。

### 验证目标

两个封闭、共位、单相单元从初始压差出发时，应满足：高压侧压力下降、低压侧压力上升、压差单调衰减；加权平均压力保持不变；首步衰减因子等于输入参数对应的解析值。

### 物理图像

基质初始处于高压、裂缝处于低压，流体通过双连续体交叉通道从基质流向裂缝。基质压力下降、裂缝压力上升；由于体系封闭，两侧最终趋向共同的 `1.5 MPa` 平衡压力。

### 为什么能够验证

算例中没有外边界流、重力、毛管力和力学耦合，唯一非零机制是交叉流。方程退化为

```text
Cm dpm/dt = -A (pm - pf)
Cf dpf/dt = +A (pm - pf)
```

因此压差解析解和守恒量都可直接写出。该设计把任何结果偏差归因到压力交叉流，而不是其他物理过程。

### 参考解和结果

输入为 [`A0_pressure_exchange.xml`](A0_pressure_exchange/A0_pressure_exchange.xml)，采用
`dt=0.02 s`、`maxTime=2 s`，运行输出归档于
`A0_pressure_exchange/runs/20260821/`。

```text
GEOS exit code                 0
time steps / cuts              100 / 0
time step                     0.02 s
numerical q                    0.8928571053
reference q                    0.8928571429
relative q error               4.21e-8
relative mean-pressure drift  8.81e-8
resolved-curve max error       4.39e-7 through t=1.16 s
last archived pressure contrast 11.97 Pa
```

采用合理时间步后，基质压力从约 `1.946 MPa` 下降、裂缝压力从约
`1.054 MPa` 上升，并在约 `1 s` 内接近 `1.5 MPa` 平衡。压力差在可分辨
区间内与后向 Euler 参考的最大相对误差为 `4.39e-7`，平均压力相对漂移为
`8.81e-8`，因此 A0 在可分辨区间内通过。A0 不验证毛管、重力或
compositional `PotGrad`。

详见 [A0 综合中文报告](A0_pressure_exchange_combined/report_zh.md)，其中分别给出等体积和不等体积子算例、推导、图片与数值结果。

## G0：重力驱替

### 验证目的

单独验证 Kazemi 型 gravity drainage pressure（GDP）的符号、大小及其加入矩阵侧势差的方向。

### 验证目标

初始两侧压力相同，但矩阵流体密度 `1100 kg/m3`、裂缝流体密度 `800 kg/m3`，重力 `g_z=-9.81 m/s2`，特征长度 `Lz=100 m`。目标平衡压差为

```text
p_m - p_f = -GDP
GDP = g_z (rho_f - rho_m) Lz / 2 = 147150 Pa
```

### 物理图像

两侧初始等压时，较重的基质流体在重力作用下向裂缝排泄，裂缝流体反向补入基质。系统建立负的 `p_m-p_f` 压差，直到重力势差与压力势差平衡。

### 为什么能够验证

两个连续体各只有一个共位单元，普通单元中心重力水头相互抵消；初始压力差也为零，所以剩余驱动力只能是 GDP。稳态时势差为零，直接得到上式。压差的符号还检查了重流体从矩阵向裂缝排泄的方向。

### 参考解和结果

输入为 [`G0_gravity_exchange.xml`](G0_gravity_exchange/G0_gravity_exchange.xml)，实际运行输出归档在 `G0_gravity_exchange/runs/20260821/`。

```text
GEOS exit code                 0
time steps / cuts              100 / 0
reference GDP                  147150.0000 Pa
final p_m - p_f               -147150.0803 Pa
absolute error                 0.0803 Pa
relative error                 5.46e-7
```

判据为相对误差 `<1e-5`、绝对误差 `<1 Pa`，均满足，因此 G0 通过。共同绝对压力偏移由两侧参考密度不同造成，不作为 GDP 误差。该案例是单相验证，不覆盖多相密度平均和 compositional Jacobian。

详见 [`G0 report`](G0_gravity_exchange/report.md)。

## I0：毛管/渗吸驱替

### 验证目的

验证毛管压力差是否真正进入 compositional 双连续体交叉流，并进一步引起组分转移和相饱和度变化。这是三个目标过程里最直接针对渗吸的案例。

### 验证目标

保持两侧压力、组分、温度、孔隙度、渗透率和相对渗透率一致，关闭重力；只让矩阵有水相 entry pressure、裂缝使用零毛管力模型。初始输出明确给出

```text
Pc,m(water) = 516014.125794 Pa
Pc,f(water) = 0 Pa
```

因此水相势差非零。目标是观察到非零、守恒一致的组分历史和相体积分数历史。

### 物理图像

基质水相具有较高毛管吸力，优先吸入水；水相从裂缝进入基质，基质气相被挤出并转移到裂缝。因此基质气相分数下降、裂缝气相分数上升，同时封闭体系中的 CO2 和 water 总量不变。

### 为什么能够验证

除了毛管压力外，其他宏观驱动力都被移除；同时直接输出 `phaseCapPressure`，可以证明毛管本构状态确实生成，而不是把初始化失败误认为物理无响应。因为没有组分源汇，矩阵和裂缝的组分变化还应等量反向。

该案例没有采用单一 ODE 的“解析瞬态解”，因为 CO2-brine flash、Brooks-Corey 毛管力和相对渗透率组合后，不能在不增加额外近似的情况下化为一个可靠的闭式解。因此这里采用“非零驱动 + 可观测动态响应 + 守恒”的过程判据。

### 实际结果

输入为 [`I0_capillary_exchange.xml`](I0_capillary_exchange/I0_capillary_exchange.xml)。正式的无外部压力边界封闭体系运行输出和日志归档在 `I0_capillary_exchange/runs/20260821/`。

```text
GEOS exit code                 0
accepted time steps / cuts     10449 / 0
matrix water capillary Pc      515979.664 -> 501057.703 Pa
fracture water capillary Pc    0 -> 0 Pa
matrix gas volume fraction     0.0896370 -> 0.0527413
fracture gas volume fraction   0.0897263 -> 0.1225718
CO2 total relative drift       -5.2352e-08
water total relative drift     +3.2311e-08
```

案例检查脚本 [`I0_capillary_exchange/scripts/check_I0.py`](I0_capillary_exchange/scripts/check_I0.py) 对上述输出独立计算通过。毛管状态、交叉流残量成对装配和接受时间步三个层次均已有证据，因此 I0 的动态渗吸过程通过。早期 `BiotPorosity` 版本因共同压力零空间秩亏而失败；正式输入改用压力依赖孔隙度物性提供封闭体系的压力储集项，不使用压力边界。

详见 [`I0 report`](I0_capillary_exchange/report.md)。此前带压力锚定的运行虽能正常退出并产生状态变化，但压力边界引入了外部质量，不能作为封闭渗吸验证；因此不替代本次封闭体系证据。

## J0：`PotGrad` 导数契约

### 验证目的

检查毛管势对矩阵和裂缝相体积分数的导数符号，尤其是矩阵项应为负、裂缝项应为正。

### 验证目标与原因

对简化势函数

```text
Phi(Sm,Sf) = pm - (a + b Sm) - pf + (a + b Sf)
```

解析导数是 `dPhi/dSm=-b`、`dPhi/dSf=+b`。中心有限差分只保留这段局部代数，排除了 flash、mobility、重力、网格和 Newton 误差，因此能直接发现导数符号或支持点归属错误。

### 物理图像

J0 不是一个运行时输运算例，而是局部势差的方向图像：增大基质相体积分数会增大基质毛管压力，从而降低交叉流势差；增大裂缝相体积分数会增大裂缝毛管压力，从而提高势差。对应曲线斜率分别为负和正。

### 结果和边界

脚本输出 `dPhi/dSm=-8000 Pa`、`dPhi/dSf=+8000 Pa`，相对误差低于 `1e-10`，J0 通过。但是它不是 GEOS kernel 的运行时测试，不能证明字段映射正确，也不能证明运行时确实走到了这段导数路径；I0 的动态历史证据独立决定渗吸过程是否通过。

详见 [`J0 report`](J0_potgrad_derivative/report.md)。

## P1：Mandel 流固耦合

### 验证目的

用 Mehrabian & Abousleiman (2014) 的双孔隙双渗透 Mandel 问题检验双连续体压力-变形耦合、Mandel-Cryer 过冲、基质/裂缝不同压力历史，以及 FIM 与 Sequential 路径的一致性。

### 物理图像

外部位移加载首先使体系近似不排水压缩，基质和裂缝压力快速升高并出现 Mandel-Cryer 过冲；随后流体扩散、储集和两连续体交换共同作用，压力逐渐排水衰减。基质与裂缝由于储集和渗透特性不同，保留不同的压力历史。

### 验证目标与原因

该问题具有文献给出的 N=2 解析解，几何、边界、物性和无量纲时间尺度在仓库中有明确来源；因此它能验证完整流固耦合响应，而不只是某个局部符号。FIM/Sequential 以及 effective/intrinsic 四个受控输入分别用于检查耦合方式和输入模式一致性。

### 当前证据和限制

本轮在 `P1_dpdp_mandel/runs/20260821/` 独立归档四个 same-source 输入，均为 GEOS exit code `0`、零切步；FIM 完成 `1247` 个时间步，Sequential 完成 `1373` 个时间步。四条主曲线均具备正确初始压力、Mandel-Cryer 过冲、平台段和后期排水趋势。

相对于 `analytical_script` 的平均/最大绝对误差如下：

```text
FIM effective   matrix   2.8445e-3 / 8.8370e-3
FIM effective   fracture 5.9141e-4 / 3.2508e-2
FIM intrinsic   matrix   2.9361e-3 / 9.0107e-3
FIM intrinsic   fracture 5.8439e-4 / 3.2589e-2
Seq effective   matrix   2.3054e-3 / 1.0457e-2
Seq effective   fracture 8.6255e-4 / 1.9231e-2
Seq intrinsic   matrix   2.3131e-3 / 1.0725e-2
Seq intrinsic   fracture 8.7378e-4 / 1.9361e-2
```

effective/intrinsic 压力历史的最大相对差小于 `2.63e-4`（FIM）和 `2.57e-4`
（Sequential）。P1 因局部压力和晚期排水误差、Sequential 分裂/时间离散误差，以及解析脚本与论文
Fig. 5c 数字化曲线的部分差异，仍属于**有条件定量验证**，不宣称为精确 release-gate 通过。

详见 [`P1 report`](P1_dpdp_mandel/report.md) 和归档的 [`Mandel workflow`](P1_dpdp_mandel/reference/DPDP_Mandel_validation_workflow.md)。文献全文/问题描述已随验证资料复制到 `P1_dpdp_mandel/problem_description/`，本轮不需要额外请求全文权限。

## 运行范围和判定原则

## 图件

每个子验证均有独立图件，归档在各自案例的 `figures/` 目录，互不混装：

- [A0 压差驱动图](A0_pressure_exchange/figures/A0_pressure_exchange.png)
- [G0 重力驱替图](G0_gravity_exchange/figures/G0_gravity_exchange.png)
- [I0 毛管/渗吸图](I0_capillary_exchange/figures/I0_capillary_exchange.png)
- [J0 `PotGrad` 导数图](J0_potgrad_derivative/figures/J0_potgrad_derivative.png)
- [P1 Mandel 对比图](P1_dpdp_mandel/figures/P1_dpdp_mandel_pressure_comparison.png)

1. A0/G0 的“通过”只对各自隔离的单一机制负责。
2. I0 的通过来自独立动态历史和守恒判据，不能由 J0 的代数通过替代；同理，J0 通过不等于运行时字段映射已被证明。
3. P1 的定量结论绑定本轮 `P1_dpdp_mandel/runs/20260821/` 输出、归档解析脚本和误差 CSV；不能用旧 `/tmp` 历史文件替代。
4. A0/G0/I0/P1 的当前运行输出分别放在各自案例的 `runs/20260821/`；不再依赖公共输出目录。
