# DPDP Mandel Problem — Fully Implicit Verification

## 物理问题

模拟 **双孔双渗（DPDP）多孔介质** 在压缩载荷下的 **Mandel 固结过程**，验证 Mandel-Cryer 效应（孔隙压力非单调响应：先升后降）。

参考论文：
> Mehrabian, A., & Abousleiman, Y. N. (2014). "Generalized Biot's theory and Mandel's problem of multiple-porosity and multiple-permeability poroelasticity." *J. Geophys. Res. Solid Earth*, 119(4), 2745–2763.

材料参数取自 GOM 页岩（Table 1）：**基质 (matrix)** + **大裂缝 (macro-fracture)** 双连续介质。

## 几何与网格

| 参数 | 值 |
|---|---|
| 尺寸 (2a × 2b) | 30 mm × 30 mm |
| 单元 | 10×1×10 C3D8 六面体 |
| 网格 | mesh1 (基质), mesh2 (裂缝, 同位) |

## 边界条件

```
     z=b ┌──────────────────┐  u_z = f(t)  (位移加载)
         │                  │  p: 零通量 (默认)
         │                  │
         │       y=0        │  u_y = 0
  x=0    │       y=b        │  u_y = 0
  u_x=0  │  ◄── 2a ───►    │  x=a
  ∂p/∂x=0│                  │  σ_xx=0 (应力自由)
         │                  │  p=0   (排水)
         │                  │
     z=0 └──────────────────┘  u_z = 0
                               p: 零通量 (默认)
```

| 面 | 位移 BC | 压力 BC |
|---|---|---|
| x=0 (xneg) | u_x = 0 (对称) | 零通量 ∂p/∂x=0 (自然 BC) |
| x=a (xpos) | 无约束 (应力自由) | p = 0 (排水) |
| y 两面 (yneg/ypos) | u_y = 0 | 零通量 (默认) |
| z=0 (zneg) | u_z = 0 | 零通量 (默认) |
| z=b (zpos) | u_z = f(t) (位移载荷) | 零通量 (默认) |

**关键说明**：
- `xneg` (x=0) 为对称面，约束 u_x 但不设置 pressure BC。GEOS 中无 pressure BC 的边界面默认为零通量 (∂p/∂x=0)，等价于不可渗透边界。
- `xpos` (x=a) 是排水边界且应力自由，仅设 p=0，**不约束位移**。

## 载荷

z=b 面施加 **随时间变化的解析位移载荷** `loadFunction`，由 Mandel 单孔隙度解析解校准 (k=8e-17, τ~6s)，位移量级 ~5×10⁻¹⁰ m。

## 材料参数 (GOM Shale, Table 1)

### 基质 (Matrix)
| 参数 | 值 | 来源 |
|---|---|---|
| 体积模量 K | 1.1 GPa | Table 1 |
| 剪切模量 G | 678.9 MPa | Table 1 |
| 颗粒体积模量 K_s | 27 GPa | Table 1 |
| 参考孔隙度 φ₀ | 0.14 | Table 1 |
| 渗透率 k | 4.93×10⁻²¹ m² (5 nd) | Table 1 |

### 裂缝 (Fracture)
| 参数 | 值 | 来源 |
|---|---|---|
| 体积模量 K | 22.5 MPa | Table 1 |
| 剪切模量 G | 13.89 MPa | Table 1 |
| 颗粒体积模量 K_s | 27 GPa | Table 1 |
| 参考孔隙度 φ₀ | 0.95 | Table 1 |
| 渗透率 k | 4.93×10⁻¹⁸ m² | ~~Table 1 (4.93e-15)~~ 调整值 |

### 流体 (水)
| 参数 | 值 |
|---|---|
| 密度 | 1000 kg/m³ |
| 粘度 | 0.001 Pa·s |
| 压缩系数 | 5.734×10⁻¹⁰ Pa⁻¹ |

### 窜流交换系数
| 参数 | 值 |
|---|---|
| Γ₁₃ | 1.67×10⁻²² Pa⁻¹·s⁻¹ |

> 注：渗透率已从论文原始值(基质 5 nd, 裂缝 5 md)调整，以平衡 Jacobian 矩阵，确保全隐式 block-by-block 耦合迭代收敛。

## GEOS 能力演示

| 能力 | 说明 |
|---|---|
| **SinglePhaseDualContinuumPoromechanics** | 双连续介质多孔弹性全耦合求解器 |
| **DualContinuumFVM** | 双连续介质流动框架，分别对基质/裂缝使用独立单相 FVM 求解器 |
| **SinglePhaseFVM** | 单相渗流有限体积法 (TPFA 离散) |
| **SolidMechanicsLagrangianFEM** | 固体力学拉格朗日有限元，准静态 (QuasiStatic) |
| **FullyImplicit coupling** | 全隐式 block-by-block 流固耦合 |
| **DualContinuumCrossFlow** | 基质-裂缝窜流交换，直接指定交换系数 Γ |
| **TimeHistory output** | 时间序列输出 (压力、位移) |
| **TableFunction** | 表格函数驱动的位移载荷 |
| **GMRES linear solver** | Krylov 子空间线性求解器 |
| **TwoPointFluxApproximation** | TPFA 通量近似 |
| **BiotPorosity** | Biot 孔隙度模型 |
| **CompressibleSinglePhaseFluid** | 可压缩单相流体模型 |

## 已知限制

1. **渗透率缩放**：论文原始 k (5×10⁻²¹ m²) 极低，Jacobian 中 K_pμ/K_pp 比值 ~24000，导致压力方程近奇异。当前裂缝 k 从 5 md 缩至 ~5×10⁻¹⁸ m² 以平衡矩阵条件数。论文精确值需 monolithic kernel 支持。
2. **顺序耦合未支持**：缺少 pressure_k 场注册，Sequential 耦合方式不可用。
3. **窜流模型简化**：当前通过 `interporosityExchangeCoefficient` 直接指定 Γ₁₃，而非用 Kazemi 形状因子公式自动计算。

## 文件列表

| 文件 | 说明 |
|---|---|
| `DPDP_Mandel_Mehrabian2014_FIM.xml` | GEOS 输入文件 |
| `README.md` | 本说明文件 |
