# TDDFT 激发态计算

含时密度泛函理论 (Time-Dependent Density Functional Theory, TDDFT) 是计算分子激发态能量和性质的主流方法。REST 程序支持基于 RI (Resolution of Identity) 加速的线性响应 TDDFT 计算，包括 TDA (Tamm-Dancoff Approximation) 和全线性响应 (Full LR) 两种方案，以及频率域的响应 TDDFT 计算。

此外，REST 还提供：

- **非限制 TDDFT (UTDDFT)**：以 `spin_polarization = true` 的非限制 (UHF/UKS) 参考态作激发态计算；
- **范围分离杂化 (RSH) 泛函**的 TDDFT（如 $\omega$B97X、CAM-B3LYP 等）；
- **SCF 波函数稳定性分析**：检测 SCF 解是否为鞍点（RHF/RKS 与 UHF/UKS 参考）；
- **激发态解析核梯度**：对指定激发态计算解析核梯度（`tddft_grad_state`），可直接用于几何优化等任务；
- **PySOC 导出**：将单重态/三重态振幅导出为 PySOC 可读的 JSON 文件，用于后续自旋-轨道耦合计算。

## 计算模式

REST 的 TDDFT 模块提供两种计算模式（输入关键字 `tddft_mode`）：

1. **本征值 TDDFT**：求解 Casida 方程，直接获得激发能 $\omega_n$ 和振子强度。支持 TDA 和 Full LR 方案。
2. **响应 TDDFT**：求解频率空间的线性方程组，获得频率依赖的动态极化率 $\alpha_{zz}(\omega)$。支持多种迭代求解器（仅 MO 实现、限制性参考）。

两种模式均基于 KS 轨道能量构建对角项，并通过 RI 加速库仑与交换矩阵-矢量积运算。XC 核 (fxc kernel) 通过 libxc 计算，支持 LDA、GGA、杂化与范围分离杂化泛函。

### 内核实现模式 `tddft_mode`

| 取值 | 说明 |
|---|---|
| `"mo"` (默认) | 预先将三中心 RI 积分变换到 MO 基，矩阵-矢量积全部在 MO 振幅空间完成。成熟、经过广泛验证的路径 |
| `"ao"` | Davidson 迭代仍在 MO 振幅空间，但每次矩阵-矢量积构造 AO 过渡密度，库仑/交换经 `ri_jk`、XC 核经数值格点批量计算。内存占用显著更低，且**三重态与 `tddft_spin = "both"` 仅在此模式支持** |

一般建议：常规单重态价激发计算用默认的 `"mo"` 模式；大体系（MO 基 RI 张量内存成为瓶颈）、三重态激发、以及 HF 参考的稳定性分析用 `"ao"` 模式。

### 非限制参考 (UTDDFT)

在 `[ctrl]` 中设置 `spin_polarization = true`（UHF/UKS 参考）即可进行非限制 TDDFT。非限制响应只有单一的自旋耦合通道（α 与 β 激发通道通过自旋无关的库仑核耦合），因此 `tddft_spin` 关键字不适用——在非限制计算中显式给出 `tddft_spin` 会直接报错。MO 与 AO 两种内核模式均支持非限制参考。

## 输入关键字

所有 TDDFT 相关的关键字均在输入卡的 `[tddft]` 区块中声明。

### 基本设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `tddft_method` | String | `"lr"` | TDDFT 方法。`"tda"` 为 Tamm-Dancoff 近似，`"lr"` 为全线性响应 |
| `tddft_mode` | String | `"mo"` | 内核实现模式。`"mo"` (MO 基 RI 张量) 或 `"ao"` (AO 过渡密度核) |
| `tddft_spin` | String | `"singlet"` | 自旋通道（仅限制性参考适用）。`"singlet"` 单重激发、`"triplet"` 三重激发（需 AO 模式）、`"both"` 单重+三重都算（需 AO 模式）。非限制参考 (`spin_polarization = true`) 下显式给出该关键字会报错 |
| `nroots` | usize | 6 | 需计算的激发态数目 |
| `tddft_use_optimized_fxc` | bool | true | 是否使用 rayon 并行的 fxc 矩阵-矢量积核 |
| `tddft_cutoff_energy` | f64 | 1.0e6 | 虚轨道能量截断 (Hartree)。能量高于此值的虚轨道不参与激发空间；非限制参考下 α/β 通道独立解析各自的截断窗口 |

### Davidson 求解器设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `davidson_tol` | f64 | 1.0e-10 | Davidson 迭代收敛阈值：残差范数 $\|r\| < \sqrt{\varepsilon}$ 且相邻迭代能量变化 $|\Delta E| < \varepsilon$ |
| `davidson_max_iter` | usize | 50 | Davidson 最大迭代次数 |
| `davidson_max_subspace` | usize | 60 | 最大子空间维度乘数 (子空间维度 = `nroots × davidson_max_subspace`)。Full LR 求解需要足够大的子空间以在首次重启前收敛，不建议调小 |

### AO 模式高级选项

以下关键字仅在 `tddft_mode = "ao"` 时生效（MO 模式下忽略）：

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `grid_batch` | bool | true | XC 核求值按格点分批执行，避免完整 AO-on-grid 张量常驻内存。以少量时间开销换取约 40% 的峰值内存下降 |
| `tddft_ao_rik_driver` | String | `"semitrans"` | 交换 K 的计算方式：`"semitrans"` 占据侧半变换收缩（精确，默认）；`"dm"` 精确批量密度驱动；`"lowrank"` 逐向量 SVD 低秩近似（阈值 `tddft_svd_tol`，一般不推荐） |
| `tddft_fxc_driver` | String | `"semitrans"` | XC 核的计算方式：`"semitrans"` 将 $C_{vir}$ 折入振幅、虚轨道侧直接与格点裸 AO 收缩（默认）；`"mo"` 缓存占据侧格点投影 + 虚轨道侧流式（MO 模式 fxc 算法）；`"dm"` 组装过渡密度的 NIMatmul 回退路径 |
| `tddft_svd_tol` | f64 | 1.0e-6 | `"lowrank"` K 驱动的相对奇异值阈值 |

### SCF 稳定性分析

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `stability` | String | `"off"` | SCF 波函数稳定性分析：`"internal"` (RHF/RKS 单重或 UHF/UKS 轨道 Hessian)、`"external"` (RHF→UHF 三重通道检查，仅 RHF/RKS)、`"full"` (internal + external)、`"auto"` (推荐，自动选择适用于当前参考的全部检查)。设为 `"off"` 以外的值时，本任务**只做稳定性分析、不做激发态计算** |
| `stability_nroots` | usize | 3 | 稳定性 Hessian 求解的最低本征值数目 |
| `stability_tol` | f64 | 1.0e-8 | 稳定性 Hessian 的 Davidson 收敛阈值 |

NOTE: 顶层 `[ctrl]` 区块中的 `check_stab` 关键字接受与 `stability` 相同的取值。两者同时存在时 `[tddft] stability` 优先。仅需稳定性分析时，输入卡无需实际准备激发态计算的其他设置。

### 激发态解析梯度

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `tddft_grad_state` | usize | 0 | 为第几个激发态计算解析核梯度 (1-based)；0 表示不计算。仅支持限制性参考；需要同一任务中先完成 TDDFT 求解（`tddft_grad_state` 不得超过 `nroots`） |

设置 `tddft_grad_state` 后，梯度任务 (`jobtype = force`) 与几何优化等任务的总梯度中会自动包含所选激发态的响应贡献。`tddft_spin = "both"` 时梯度取单重态通道。

### PySOC 导出

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `pysoc` | bool | false | 将 TDDFT 结果导出到 `rest_pysoc_export.json`（PySOC 格式），用于自旋-轨道耦合后续计算。要求限制性参考且 `tddft_spin = "both"`（自旋-轨道耦合同时需要单重与三重振幅） |

### 响应 TDDFT 设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `response_tddft` | bool | false | 是否启用响应 TDDFT 计算（仅 MO 实现、限制性参考） |
| `response_tddft_solver` | String | `"klopper"` | 响应方程求解器。`"pople"`, `"gmres"`, `"klopper"`, `"dense"` 四选一 |
| `response_tddft_tol` | f64 | 1.0e-6 | 迭代求解器收敛阈值 |
| `response_tddft_max_iter` | usize | 200 | 迭代求解器最大迭代次数 |
| `external_field_freq` | f64 | 0.5 | 外场频率 $\omega$ (Hartree) |
| `lifetime_gamma` | f64 | 0.001 | 寿命展宽 $\gamma$ (Hartree) |

### 响应 TDDFT 空间格点采样

当启用 `response_tddft = true` 时，以下关键字用于生成空间采样格点以导出极化密度：

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `response_tddft_x_start` | f64 | 0.0 | X 方向起始坐标 (Bohr) |
| `response_tddft_x_end` | f64 | 1.0 | X 方向终止坐标 (Bohr) |
| `response_tddft_x_points` | usize | 2 | X 方向采样点数 |
| `response_tddft_y_start` | f64 | 0.0 | Y 方向起始坐标 (Bohr) |
| `response_tddft_y_end` | f64 | 1.0 | Y 方向终止坐标 (Bohr) |
| `response_tddft_y_points` | usize | 2 | Y 方向采样点数 |
| `response_tddft_z_start` | f64 | 0.0 | Z 方向起始坐标 (Bohr) |
| `response_tddft_z_end` | f64 | 1.0 | Z 方向终止坐标 (Bohr) |
| `response_tddft_z_points` | usize | 2 | Z 方向采样点数 |

NOTE: 采样格点按外循环 X、中循环 Y、内循环 Z 的顺序生成。当某方向点数 ≤1 时不实际采样。

### FEAST 求解器设置

当 `tddft_feast_solver = true` 时，以下关键字生效（围道积分求解器仅 MO 模式、限制性参考支持）：

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `tddft_feast_solver` | bool | false | 是否使用 FEAST 围道积分求解器 |
| `tddft_feast_eigenrange_min` | f64 | 0.0 | 能量窗口下界 (Hartree) |
| `tddft_feast_eigenrange_max` | f64 | 0.5 | 能量窗口上界 (Hartree) |
| `tddft_feast_m_expected` | usize | 20 | 窗口内预期本征值数目 |
| `tddft_feast_max_iter` | usize | 30 | FEAST 外迭代最大次数 |
| `tddft_feast_tol` | f64 | 1.0e-8 | FEAST 收敛阈值 |
| `tddft_feast_gmres_restart` | usize | 200 | GMRES 重启维度 |
| `tddft_feast_gmres_max_iter` | usize | 500 | GMRES 内迭代最大次数 |
| `tddft_feast_cg_max_iter` | usize | 100 | CG 内迭代最大次数 |
| `tddft_feast_cg_tol` | f64 | 1.0e-8 | CG 收敛阈值 |
| `tddft_feast_init_guess_type` | String | `"random"` | 初始猜测类型 |
| `tddft_feast_gaussian_width_factor` | f64 | 0.5 | 高斯宽度因子 |

## 用法示例

### 示例 1：H2 分子的 TDA 单重激发

```toml
[ctrl]
     print_level =               2
     num_threads =               4
     xc =                        "pbe0"
     basis_path =                "def2-tzvp"
     auxbas_path =               "def2-tzvp-rifit"
     eri_type =                  "ri-v"
     use_ri_symm =               true
     charge =                    0.0
     spin =                      1.0
     spin_polarization =         false
     mixer =                     "diis"
     max_scf_cycle =             100
     scf_acc_rho =               1.0e-8
     initial_guess =             "sad"

[tddft]
     tddft_method =              "tda"
     nroots =                    3

[geom]
     name = "H2"
     unit = "angstrom"
     position = '''
        H  0.0  0.0  0.0
        H  0.0  0.0  0.74
     '''
```

程序首先完成 SCF 收敛，随后调用 TDDFT 模块计算最低 3 个激发态的能量和振子强度。输出将逐态打印激发能 (eV)、振子强度、跃迁偶极矩平方等信息，并把激发能/振子强度汇总写入 `rest_results.json` 的 `"tddft"` 字段。

### 示例 2：三重激发（AO 模式）

```toml
[tddft]
     tddft_method =              "tda"
     tddft_mode =                "ao"
     tddft_spin =                "triplet"
     nroots =                    6
```

### 示例 3：SCF 稳定性分析

```toml
[tddft]
     stability =                 "auto"
     stability_nroots =          3

[ctrl]
     check_stab =                "auto"
```

NOTE: `stability` 与 `check_stab` 任一非 `"off"` 即触发稳定性分析（前者优先）；本任务只做稳定性检查，不计算激发态。最低本征值 $\lambda_{\min} < -10^{-5}$ 即报告不稳定，结果同时打印到输出并写入 `rest_results.json` 的 `"stability"` 字段。

### 示例 4：激发态解析梯度

```toml
[tddft]
     tddft_method =              "tda"
     nroots =                    3
     tddft_grad_state =          1
```

在梯度任务 (`jobtype = force`) 中，程序先求解 TDDFT，再将第 1 个激发态的响应贡献加入基态梯度，输出总梯度。该梯度同样可用于几何优化。

### 示例 5：非限制 TDDFT

```toml
[ctrl]
     spin_polarization =         true
     spin =                      1.0

[tddft]
     tddft_method =              "tda"
     nroots =                    6
```

非限制计算输出 α/β 两个自旋通道的激发（如 `#3a->#5b` 表示 α 轨道 3 到 β 轨道 5 的激发），激发能/振子强度的后处理约定与 PySCF UHF-TDDFT 一致。

## 注意事项

- **RI 加速前提**：TDDFT 模块依赖 RI 积分加速。输入卡中必须提供 `auxbas_path` 并设置 `eri_type = "ri-v"`。RSH 泛函的短程交换积分在 SCF 阶段自动构建，无需额外设置。
- **求解器选择**：对小型体系 (限制性参考激发空间维度 ≤15；非限制 MO 模式 TDA ≤15、Full LR ≤80)，程序自动使用稠密对角化 (Full LR 为非厄米 $[\mathbf{A}\ \mathbf{B}; -\mathbf{B}\ -\mathbf{A}]$ 直接对角化)；中等以上体系默认使用 Davidson 迭代求解器 (AO 模式为批量接口)；特定能量窗口计算可启用 FEAST 求解器（仅限制性 + MO 模式）。
- **单重/三重激发**：限制性参考通过 `tddft_spin` 控制；三重态与 `"both"` 需要 `tddft_mode = "ao"`。非限制参考只有单一自旋耦合通道，不接受 `tddft_spin`。
- **稳定性分析**：与激发态计算互斥；需要 DFT 参考提供数值格点（HF 参考无需格点，自动只算 RI J/K 部分）；ROHF 参考不支持；实→复稳定性与 UHF→GHF 外稳定性尚未实现。
- **解析梯度**：仅限制性参考；需先在同行任务中完成 TDDFT 求解；梯度为 RI 近似下的解析梯度。
- **响应 TDDFT**：响应模式的线性系统为 4 分量非厄米方程组，维度为本征值模式的 4 倍。Klopper 子空间求解器是推荐选择。
