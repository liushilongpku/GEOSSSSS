# CO2-Brine WAG 双重介质流固耦合范例

## 范例定位

本目录是从远程工作站案例整理出的输入文件范例。它保留了运行该案例所需的
`input.xml`、网格、流体参数文件和相对渗透率/初始化表格，但**不包含运行输出**，例如
`output/`、`run.log`、`run.pid`、`run.exit`、`run_metadata.txt` 或结果图片。

当前本地范例目录为：

```text
/home/lsl/codes/GEOSSSSS/inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkHMCompositional/co2_brine_wag/
```

主输入文件为：

```text
co2_brine_wag/input.xml
```

## 模拟问题

该输入文件模拟四层储层/盖层/基底/断层结构中的 **CO2-盐水两相 WAG（Water Alternating Gas）注入与生产过程**。
WAG 的含义是向储层交替注入 CO2 和水，通过周期性改变注入相来研究 CO2 封存、盐水驱替、压力传播、
相态变化以及基质-裂缝之间的物质交换。

案例包含两个重合的网格体系：

- `mesh1`：基质连续介质，包含 aquifer、cap、reservoir、basement 和 fault 等区域；
- `mesh2`：裂缝连续介质，主要对 reservoir 区域进行显式裂缝流动计算。

注入井位于裂缝储层区域 `reservoirF`，生产井位于基质储层区域 `reservoirM` 和裂缝储层区域
`reservoirF`。模型因此可以同时观察：

- 基质和裂缝中的压力、CO2 饱和度及水饱和度；
- CO2 在盐水中的溶解和盐水中的水蒸气组分；
- CO2 注入、水注入和生产边界的周期变化；
- 基质-裂缝之间的跨连续介质传输；
- 重力驱替、毛管压力和滞回相对渗透率对 WAG 响应的影响；
- 固体变形、应力变化和孔隙度/渗透率随流固耦合状态的变化。

该案例不是单一均质块体的解析验证算例，而是面向多区域、多相态、双重介质和长期 WAG 响应的应用型范例。

## 主要物理模型

`input.xml` 中使用的主要求解器和本构模型包括：

- `CompositionalMultiphaseFVM`：基质和裂缝中的 CO2-盐水组分多相流；
- `CompositionalMultiPhaseDualContinuumFVM`：基质-裂缝双连续介质流动及交换；
- `SolidMechanicsLagrangianFEM`：准静态固体力学；
- `CompositionalMultiphaseDualContinuumPoromechanics`：多相流-固体力学耦合；
- `CO2BrinePhillipsFluid`：CO2-盐水相态与 PVT 模型；
- `SimpleGravityDrainagePressure`：重力驱替压力模型；
- `TableRelativePermeabilityHysteresis`：基于表格的 Killough 滞回相对渗透率；
- `BrooksCoreyCapillaryPressure`：Brooks-Corey 毛管压力；
- `ElasticIsotropic`、`BiotPorosity`、`ConstantPermeability` 和 `PorousElasticIsotropic`：多区域孔隙弹性材料；
- `DualContinuumCrossFlow`：基质-裂缝跨连续介质交换。

## WAG 工况

该范例属于 `5.15hyst_base_stopprod` 基准案例，具有以下工况特征：

- WAG 单段周期：3 个月；
- CO2 注入通量：`-20.0`；
- 水注入通量：`-48.89`；
- 水/CO2 质量比：约为 1；
- 注入边界和生产压力边界持续到 `315360000 s`，即约 10 年；
- 未显式设置 `newtonMinIter`，继承当前 GEOS 程序默认值 `1`；这与远程 `20260807_miniter1` 版本的行为一致；
- 双连续介质裂缝体积分数：由 XML 中的 `fractureVolumeFraction` 设置，当前输入为 `0.1`；
- 开启重力驱替和滞回相对渗透率。

具体数值、分段时间函数、区域材料参数和边界集合以 `input.xml` 为准。README 不重复抄录全部 XML 参数，
避免说明文字与输入文件发生漂移。

## 输入文件结构

```text
co2_brine_wag/
├── input.xml
├── reservoir_4layer_400x200_fault_wells.vtk
└── tables/
    ├── co2flash.txt
    ├── pvtgas.txt
    ├── pvtliquid.txt
    ├── fluid_phaseModel2_SpanWagnerCO2Density_table.csv
    ├── pressure.csv
    ├── sigma_h.csv
    ├── sigma_v.csv
    ├── x.csv
    ├── y.csv
    ├── z.csv
    ├── drainage_Sw.txt
    ├── drainage_Sg.txt
    ├── drainage_krw.txt
    ├── drainage_krg.txt
    ├── imbibition_Sw.txt
    ├── imbibition_Sg.txt
    ├── imbibition_krw.txt
    └── imbibition_krg.txt
```

表格已经集中放入 `tables/`，并且 `input.xml` 中的引用已经按本地目录布局调整为
`tables/<filename>`。因此，运行时应从本目录或使用绝对路径指定 XML，不能把 XML 单独复制到其他位置后
再期待相对路径仍然有效。

## 来源与原始文件位置

该范例来自 `dualporo_HM` 项目的 WAG 研究案例。远程工作站上的原始案例目录为：

```text
/data/datafile/lsl/project/dualporo_HM/fracture_no_hyst_series_20yr_20260807_miniter1/5.15hyst_base_stopprod/
```

原始主输入文件为：

```text
/data/datafile/lsl/project/dualporo_HM/fracture_no_hyst_series_20yr_20260807_miniter1/5.15hyst_base_stopprod/input.xml
```

远程工作站 SSH 别名为 `server1`，主机名为 `hello-Precision-7875-Tower`。本地 WSL 中可用以下命令查询原始文件：

```bash
ssh server1 'ls -la /data/datafile/lsl/project/dualporo_HM/fracture_no_hyst_series_20yr_20260807_miniter1/5.15hyst_base_stopprod/'
ssh server1 'sed -n "1,220p" /data/datafile/lsl/project/dualporo_HM/fracture_no_hyst_series_20yr_20260807_miniter1/5.15hyst_base_stopprod/input.xml'
```

同名的正式基准运行系列还存在于：

```text
/data/datafile/lsl/project/dualporo_HM/fracture_no_hyst_series_20yr_20260719/5.15hyst_base_stopprod/
```

两套 `input.xml` 的主要差异是 Newton 非线性迭代参数：

- `20260719` 正式版本在三个流动/双连续介质非线性参数块中显式使用 `newtonMinIter="0"`；
- `20260807_miniter1` 版本在对应位置显式使用 `newtonMinIter="1"`；
- 本地 `co2_brine_wag/input.xml` 未显式写该属性，因此使用程序默认值 `1`，无需重复添加。

本地 `co2_brine_wag/input.xml` 取自后者。远程原始案例的运行输出仍保留在远程目录的 `output/` 中，
但本地范例有意不复制这些输出文件。

项目中的案例与数据索引可进一步查询：

```text
/home/lsl/project/wagInFracturedRock/case_and_data_locations.md
```

其中记录了 WAG 频率系列、水/CO2 质量比系列、远程原始路径以及当前论文使用的案例关系。

## 运行前注意事项

输入包已通过当前本地 GEOS 的输入加载验证，但加载验证不等于完整长期计算已经在本地完成。运行前应确认：

- 使用与该输入兼容的 GEOS 构建；
- 从本目录运行，或保证 XML 的相对路径基准目录正确；
- 输出目录放在本目录之外或另建独立目录，避免污染输入范例；
- 长期计算需要较大的内存和磁盘空间；
- 当前目录只保留输入文件，不把远程 `output/`、日志或运行状态文件当作输入的一部分。

示例加载验证命令：

```bash
build/bin/geosx -v \
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkHMCompositional/co2_brine_wag/input.xml \
  -o /tmp/geos_co2_brine_wag_validation
```

完整运行时应使用新的输出目录，例如：

```bash
build/bin/geosx \
  -i inputFiles/AsmokeTestForModifiedByL/dpdk/dpdkHMCompositional/co2_brine_wag/input.xml \
  -o /tmp/geos_co2_brine_wag_run
```
