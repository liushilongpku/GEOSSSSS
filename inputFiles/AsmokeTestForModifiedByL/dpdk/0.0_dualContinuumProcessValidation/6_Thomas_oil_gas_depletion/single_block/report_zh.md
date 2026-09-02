<!-- 目的：记录 Thomas 1983 case 6 气油单块模型的原文映射、可复现结果和剩余限制。 -->

# Thomas 1983 Case 6 单块气油重力排驱

## 1. 验证对象

本模型以一个基质单元和一个共置裂缝单元表示 Thomas（1983）的 10 ft 气油块，不通过网格细分、
形状因子、相渗端点或目标采收率反标进行拟合。

| 项目 | 设置 | 依据 |
| --- | --- | --- |
| 网格 | 基质 `1x1x1`，裂缝 `1x1x1` | Thomas 单胞模型 |
| 块尺寸 | 3.05 m 立方体 | Thomas 10 ft 块 |
| 基质孔隙度/渗透率 | 0.30 / 1 md | Thomas Table 1 |
| 初始压力 | 38.3327805 MPa | Table 2/4 的 5545 psig 绝对压力 |
| 裂缝衰竭 | 0.75 psi/day，2.5 年 | Thomas case 6 正文 |
| 形状因子 | 0.215278208 m^-2 | 原值 `0.02 ft^-2` 的 SI 换算 |
| 相渗 | GEOS 表格 modified Stone-II，气油支路使用初压 Table 4 `krg/kro` | 复用现有三相相渗模型 |
| 伪毛管压力 | 二维 `Pcgo(Sg,P)` | Table 3、Eq. (28) 与 10-ft 垂向平衡重建 |
| 显式 GDP | canonical GDP-on 使用半块高 `Lz/2` | Thomas Eq. (38) |

Thomas 正文一处写 5540 psig，而 Table 2/4 使用 5545 psig。模型采用 5545 psig，以与 PVT 和伪函数一致。

## 2. 统计口径

Thomas 按地面油计量。GEOS 的对应量是矩阵油组分质量采收率：

```text
R_o = 1 - matrix_oil_component_mass(t) / matrix_oil_component_mass(0)
```

该库存与 Thomas 的 `phi*b_o*S_o = phi*S_o/B_o` 成正比。`1-So/So0` 只是饱和度诊断；canonical
GDP-on 在 2.5 年的该诊断值为 `49.4061%`，不能代替地面油采收率 `45.5423%`。精细垂向旁证算例
按质量口径给出 `46.0030%`，也证明 Thomas 的 46% 不是简单饱和度降幅。

## 3. 压力相关垂向平衡伪毛管压力

Thomas Eq. (28) 缩放的是 Table 3 局部岩石毛管压力：

```text
Pcgo_local(Sg,P) = sigma(P)/sigma_I * Pcgo_local(Sg,P_I)
```

Table 4 已经是包含 10-ft 垂向重力积分的块平均伪函数。将整条 Table 4 直接乘
`sigma(P)/sigma_I` 会把其中的重力偏移也错误缩放；完全不考虑压力变化同样不正确。

`scripts/generate_vertical_equilibrium_pseudocapillary.py` 在每个 Table 2 压力执行以下步骤：

1. 用 Eq. (28) 缩放 Table 3 局部 `Pcgo(Sg)`。
2. 用当前 `Bo/Bg/Rs` 计算油气密度。
3. 在完整 10-ft 块高上重新积分垂向平衡。
4. 输出普通二维 `TableFunction` 所需的 `Sg` 轴、绝对压力轴和 `Pcgo` 值。

独立积分在初压对 Table 4 六个内部点的最大误差为 `0.006699 psi`。初压 Table 4 负责锚定参考曲线，
压力变化完全来自上述物理重建，不使用采收率拟合。

`PressureScaledTableCapillaryPressure` 新增 `pressureDependentTableName`，双线性计算 `Pcgo(Sg,P)` 和
`dPcgo/dSg`；三相模式仍保留水油毛管压力。该选项与旧的 `pressureScalingTableName` 互斥；未配置二维表时
旧标量路径保持不变。

## 4. 交换闭合与 GDP

方向相关 PPU 交换规则为：

- 基质到裂缝使用当前基质 mobility 和组成。
- 裂缝到基质使用真实裂缝 PVT、组成和相覆盖率。
- 反向有效相渗按相序 `{ oil, gas, water }` 设置为 `{ -1, 0.42, 0.03 }`；负油值表示使用当前基质 `kro`。
- 纯气裂缝的油水覆盖率为零，因此不会构造非物理油水回灌。

Thomas GDP 只加到油相交换势：

```text
GDP_o = (rho_o-rho_g)*abs(g_z)*Lz/2
GDP_g = GDP_w = 0
```

二维 Table 4 伪闭合与半高 GDP 在当前 GEOS 映射中都需要。关闭 GDP 明显欠排驱；使用完整块高则过驱到
`54.6939%`，因此不能用全高修正终值。

## 5. Canonical 结果

| 时间 / 年 | GDP-on / % | GDP-off / % | Thomas Fig. 4 / % | GDP-on 减 Thomas / pp |
| ---: | ---: | ---: | ---: | ---: |
| 0.25 | 17.4678 | 8.5879 | 12.50 | +4.9678 |
| 0.50 | 27.3364 | 13.2576 | 27.00 | +0.3364 |
| 0.75 | 33.5632 | 16.0687 | 35.25 | -1.6868 |
| 1.00 | 37.7124 | 17.6936 | 40.00 | -2.2876 |
| 1.50 | 42.5042 | 18.7350 | 44.15 | -1.6458 |
| 2.00 | 44.7305 | 18.7370 | 46.00 | -1.2695 |
| 2.50 | **45.5423** | **18.7370** | **46.00** | **-0.4577** |

GDP-on 终态 `So/Sg/Sw=0.40507/0.39479/0.20014`，GDP-off 为
`0.60449/0.19537/0.20014`。两次运行均为 3,963 步、零切步、零丢弃 nonlinear iteration。

0.5 年和 2.5 年分别与 Thomas 相差 `+0.3364 pp` 和 `-0.4577 pp`，落在 Fig. 4 约 `+/-1 pp` 的读图
不确定性内。但 0.25 年高约 `4.95 pp`，0.75 至 2.0 年低约 `1.29--2.31 pp`，因此不能声称整条瞬态曲线
严格重合。

## 6. 独立 Oracle

`scripts/thomas_single_cell_oracle.py` 独立积分黑油库存和单胞交换方程，并使用同一物理定义生成的二维 VE 曲线及
半高 GDP。裂缝 `Pc=0` 时 2.5 年结果为 `45.4514%`，与 GEOS 相差 `0.0909 pp`。11 个季度采样点上的
最大绝对差同样为 `0.0909 pp`，说明二者在完整历史上闭合，而不仅是终点偶合。

裂缝取端点伪毛管压力的诊断结果为 `65.3323%`，不是 canonical 边界条件。

## 7. 复现命令

在本目录运行：

```bash
python3 scripts/generate_vertical_equilibrium_pseudocapillary.py
python3 scripts/prepare_thomas_single_block.py
/home/lsl/codes/GEOSSSSS/build/bin/geosx \
  -i thomas_singleblock_gas_oil_gravity_drainage_gdp_on.xml \
  -o /tmp/thomas_singleblock_canonical_gdp_on
python3 scripts/analyze_results.py \
  --output /tmp/thomas_singleblock_canonical_gdp_on \
  --analysis-dir /tmp/thomas_analysis_canonical_gdp_on \
  --figures-dir /tmp/thomas_figures_canonical_gdp_on
python3 scripts/thomas_single_cell_oracle.py --output analysis/thomas_single_cell_oracle_ve.csv
```

所有运行输出均写入独立 `/tmp` 目录。

## 8. 已排除解释与剩余限制

- `1x1x2/4/8` 的 2.5 年结果为 `48.4454%/50.5551%/51.0888%`，不是网格收敛；每个子单元重复完整块尺度闭合。
- 三倍形状因子在 0.5 年已达 `44.0406%`，早期严重过驱，不能用于拟合 46%。
- 解除裂缝持续组分约束、改变初始 `Pcgo` 或改用气压主变量均不能解释旧结果差异。
- 交换 transmissibility 仅由基质渗透率决定；Jacobian 为 `dT/dP_matrix`，`dT/dP_fracture=0`。
- 三相相渗复用 GEOS modified Stone-II；它沿本算例近束缚水轨迹逼近 Table 4 气油伪函数，不引入专用相渗源码。
- Fig. 4 中间点来自图像读取；0.25 年和中期残余瞬态差异目前尚无原文证据可进一步消除。
