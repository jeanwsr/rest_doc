# TDDFT 激发态计算

含时密度泛函理论 (Time-Dependent Density Functional Theory, TDDFT) 是计算分子激发态能量和性质的主流方法。REST 程序支持基于 RI (Resolution of Identity) 加速的线性响应 TDDFT 计算，包括 TDA (Tamm-Dancoff Approximation) 和全线性响应 (Full LR) 两种方案，以及频率域的响应 TDDFT 计算。

## 计算模式

REST 的 TDDFT 模块提供两种计算模式：

1. **本征值 TDDFT**：求解 Casida 方程，直接获得激发能 $\omega_n$ 和振子强度。支持 TDA 和 Full LR 方案。
2. **响应 TDDFT**：求解频率空间的线性方程组，获得频率依赖的动态极化率 $\alpha_{zz}(\omega)$。支持多种迭代求解器。

两种模式均基于 KS 轨道能量构建对角项，并通过 RI 加速库仑与交换矩阵-矢量积运算。XC 核 (fxc kernel) 通过 libxc 计算，支持 LDA、GGA 和杂化泛函。

## 输入关键字

所有 TDDFT 相关的关键字均在输入卡的 `[tddft]` 区块中声明。

### 基本设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `tddft_method` | String | `"lr"` | TDDFT 方法。`"tda"` 为 Tamm-Dancoff 近似，`"lr"` 为全线性响应 |
| `tddft_spin` | String | `"singlet"` | 自旋多重度。`"singlet"` 为单重激发，`"triplet"` 为三重激发 |
| `nroots` | usize | 6 | 需计算的激发态数目 |
| `tddft_use_optimized_fxc` | bool | true | 是否使用 rayon 并行的 fxc 矩阵-矢量积核 |
| `tddft_cutoff_energy` | f64 | 1.0e6 | 虚轨道能量截断 (Hartree)。能量高于此值的虚轨道不参与激发空间 |

### Davidson 求解器设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `davidson_tol` | f64 | 1.0e-6 | Davidson 迭代收敛阈值 |
| `davidson_max_iter` | usize | 50 | Davidson 最大迭代次数 |
| `davidson_max_subspace` | usize | 8 | 最大子空间维度乘数 (子空间维度 = `nroots × davidson_max_subspace`) |

### 响应 TDDFT 设置

| 关键字 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `response_tddft` | bool | false | 是否启用响应 TDDFT 计算 |
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

当 `tddft_feast_solver = true` 时，以下关键字生效：

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

以下以一个 H2 分子的 TDA 单重激发计算为例：

```toml
[ctrl]
     print_level =               2
     num_threads =               4
     xc =                        "pbe"
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
     tddft_spin =                "singlet"
     nroots =                    3
     davidson_tol =              1.0e-6
     davidson_max_iter =         50

[geom]
     name = "H2"
     unit = "angstrom"
     position = '''
        H  0.0  0.0  0.0
        H  0.0  0.0  0.74
     '''
```

程序首先完成 SCF 收敛，随后调用 TDDFT 模块计算最低 3 个激发态的能量和振子强度。输出将逐态打印激发能 (eV)、振子强度、跃迁偶极矩平方等信息。

若需进行全线性响应 (Full LR) 计算，将 `tddft_method` 改为 `"lr"` 即可。响应 TDDFT 模式下，设置 `response_tddft = true` 并根据需要指定 `external_field_freq` 和求解器类型。

## 注意事项

- **RI 加速前提**：TDDFT 模块依赖 RI 积分加速。输入卡中必须提供 `auxbas_path` 并设置 `eri_type = "ri-v"`。
- **求解器选择**：对小型体系 (激发空间维度 ≤15)，程序自动使用稠密对角化；对中等以上体系，默认使用 Davidson 迭代求解器；对于特定的能量窗口计算，可启用 FEAST 求解器。
- **单重/三重激发**：通过 `tddft_spin` 控制。
- **响应 TDDFT**：响应模式的线性系统为 4 分量非厄米方程组，维度为本征值模式的 4 倍。Klopper 子空间求解器是推荐选择。
