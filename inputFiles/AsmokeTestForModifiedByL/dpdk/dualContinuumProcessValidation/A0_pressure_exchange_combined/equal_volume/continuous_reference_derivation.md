<!-- Purpose: derive the A0 continuous-time pressure-exchange reference from the reduced two-storage model. -->

# A0 连续时间基准推导

本文档推导 A0 压差驱动算例使用的连续时间参考解。该解来自简化的两个封闭储集体压力交换方程，用于描述基质和裂缝之间仅由压力差驱动的交叉流。

## 1. 简化模型

设基质和裂缝是两个封闭储集体，仅通过交叉流交换流体：

$$
C_m\frac{dp_m}{dt}=-A(p_m-p_f)
$$

$$
C_f\frac{dp_f}{dt}=+A(p_m-p_f)
$$

其中：

- $p_m,p_f$：基质和裂缝压力；
- $C_m,C_f$：两侧压力储集系数；
- $A$：交叉流传递系数；
- $A(p_m-p_f)$：从基质流向裂缝的交换量。

这里假设没有外部边界流、重力、毛管压力和力学耦合，并将储集系数和交叉流系数视为常数。

## 2. 储集量守恒

将两条方程相加：

$$
C_m\frac{dp_m}{dt}+C_f\frac{dp_f}{dt}=0
$$

因此：

$$
C_mp_m+C_fp_f=\text{constant}
$$

定义储集加权平均压力：

$$
\bar p=\frac{C_mp_m+C_fp_f}{C_m+C_f}
$$

则 $\bar p$ 随时间保持不变。

## 3. 压差方程

定义基质与裂缝的压力差：

$$
\Delta p=p_m-p_f
$$

对其求导并代入原方程：

$$
\begin{aligned}
\frac{d\Delta p}{dt}
&=\frac{dp_m}{dt}-\frac{dp_f}{dt}\\
&=-\frac{A}{C_m}(p_m-p_f)-\frac{A}{C_f}(p_m-p_f)\\
&=-A\left(\frac{1}{C_m}+\frac{1}{C_f}\right)\Delta p
\end{aligned}
$$

定义衰减率：

$$
\lambda=A\left(\frac{1}{C_m}+\frac{1}{C_f}\right)
$$

于是得到：

$$
\frac{d\Delta p}{dt}=-\lambda\Delta p
$$

## 4. 连续时间解

分离变量并积分：

$$
\frac{d\Delta p}{\Delta p}=-\lambda\,dt
$$

$$
\ln\frac{\Delta p(t)}{\Delta p_0}=-\lambda t
$$

因此：

$$
\boxed{\Delta p(t)=\Delta p_0e^{-\lambda t}}
$$

其中：

$$
\Delta p_0=p_{m,0}-p_{f,0}
$$

结合平均压力守恒和压差定义，可得两侧压力：

$$
\boxed{
p_m(t)=\bar p+\frac{C_f}{C_m+C_f}\Delta p_0e^{-\lambda t}
}
$$

$$
\boxed{
p_f(t)=\bar p-\frac{C_m}{C_m+C_f}\Delta p_0e^{-\lambda t}
}
$$

## 5. A0 的参数代入

A0 中两侧储集系数相同：

$$
C_m=C_f=C
$$

初始压力为：

$$
p_{m,0}=2.0\ \mathrm{MPa},\qquad p_{f,0}=1.0\ \mathrm{MPa}
$$

因此：

$$
\bar p=\frac{p_{m,0}+p_{f,0}}{2}=1.5\ \mathrm{MPa}
$$

$$
\Delta p_0=p_{m,0}-p_{f,0}=1.0\ \mathrm{MPa}
$$

A0 使用的简化系数为：

$$
A=1.2\times10^{-9},\qquad C=4.0\times10^{-10}
$$

所以：

$$
\lambda=\frac{2A}{C}
=\frac{2(1.2\times10^{-9})}{4.0\times10^{-10}}
=6\ \mathrm{s^{-1}}
$$

最终得到图中使用的连续时间参考：

$$
\boxed{p_m(t)=1.5+0.5e^{-6t}\ \mathrm{MPa}}
$$

$$
\boxed{p_f(t)=1.5-0.5e^{-6t}\ \mathrm{MPa}}
$$

## 6. 物理含义

初始时基质压力高于裂缝压力，因此流体从基质流向裂缝。基质压力下降，裂缝压力上升，压力差按指数规律衰减。由于体系封闭，两侧最终趋近共同平衡压力：

$$
p_m=p_f=1.5\ \mathrm{MPa}
$$

该连续解是简化线性模型的数学参考，不是从 GEOS 源码中直接提取的公式，也不是文献中的外部基准解。

## 7. 与 GEOS 结果的关系

GEOS 的 A0 结果还包含实际离散和物性更新的影响。因此该参考主要用于检查：

1. 压差衰减方向是否正确；
2. 基质与裂缝压力变化是否方向相反；
3. 两侧是否趋近共同平衡压力；
4. 在时间步足够小时，GEOS 曲线是否接近连续时间解。

A0 图同时绘制向后欧拉参考，因为 GEOS 使用有限时间步推进；向后欧拉参考用于区分连续模型误差和时间离散误差。
