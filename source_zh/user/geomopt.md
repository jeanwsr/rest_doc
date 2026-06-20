# 构型优化

构型优化 (Geometry Optimization) 是在势能面上搜索稳定点（能量极小点或过渡态鞍点）的过程。REST 程序通过迭代地计算分子受力（能量梯度）并据此调整原子坐标，逐步逼近目标构型，直至满足预设的收敛标准。

REST 支持两种构型优化引擎：

1. **geomeTRIC**（缺省，推荐）：基于 Python 的 [geomeTRIC](https://geometric.readthedocs.io/en/latest/) 优化库，通过 PyO3 接口与 REST 集成。支持多种内坐标系统、过渡态搜索、IRC 路径计算、约束优化（固定原子）和频率分析。
2. **L-BFGS**：基于 Rust 的 L-BFGS 算法，仅使用数值梯度，不支持约束优化和过渡态搜索，适用于简单的测试场景。

## 基本设置

通过 `[ctrl]` 区块中的关键词设置构型优化任务：

- `job_type`：设为 `"opt"` 以启动构型优化。等价设置有 `"geometry optimization"`、`"relax"`、`"geom_opt"` 等。
- `opt_engine`：优化引擎选择。缺省为 `"geometric_pyo3"`（geomeTRIC 引擎），也可设为 `"lbfgs"`。
- `numerical_force`：是否使用数值力（有限差分）代替解析梯度。对于不具备解析梯度的后自洽场方法（如 XYG3、RPA 等），需设为 `true`。缺省为 `false`。

一个基本的构型优化输入示例如下：

```toml
[ctrl]
    job_type =                  "opt"
    xc =                        "x3lyp"
    empirical_dispersion =      "d3bj"
    basis_path =                "cc-pVDZ"
    auxbas_path =               "def2-universal-jkfit"
    charge =                    0.0
    spin =                      1.0
    spin_polarization =         false

[geom]
    name = "CO"
    unit = "angstrom"
    position = """
       C  0.00000000000      0.0000000000      0.0000000000
       O  1.20000000000      0.0000000000      0.0000000000
    """
```

对于使用双杂化泛函等没有解析梯度的方法，需要开启数值力：

```toml
[ctrl]
    xc =                        "xyg3"
    job_type =                  "opt"
    numerical_force =           true
    basis_path =                "cc-pVDZ"
    auxbas_path =               "def2-universal-jkfit"
```

## 固定原子（约束优化）

在构型优化中，用户可以固定部分原子的坐标，仅优化其余原子的位置。这在表面吸附、QM/MM 等计算场景中非常实用。

固定原子通过在 `[geom]` 区块的 `position` 中添加标记实现，格式为：

```text
<元素> <固定标记> <x> <y> <z>
```

其中，固定标记为 `0` 表示该原子**固定**（不参与优化），`1` 表示该原子**自由**（参与优化）。如果不提供固定标记（即标准的四列格式），则所有原子默认为自由原子。

```{note}
固定原子功能仅在 geomeTRIC 引擎 (`opt_engine = "geometric_pyo3"`) 下可用。L-BFGS 引擎不支持约束优化。
```

## geomeTRIC 优化器参数

geomeTRIC 引擎的参数通过输入卡中的 `[ctrl.geometric_pyo3]` 子区块声明。对于大多数计算，使用缺省参数即可，无需额外设置此区块。

### 收敛标准

| 关键词 | 类型 | 缺省值 | 单位 | 说明 |
|---|---|---|---|---|
| `maxiter` | i32 | 300 | — | 最大优化步数 |
| `convergence_energy` | f64 | 1.0e-6 | Hartree | 相邻两步能量变化的收敛阈值 |
| `convergence_grms` | f64 | 3.0e-4 | Hartree/Bohr | 梯度 RMS 的收敛阈值 |
| `convergence_gmax` | f64 | 4.5e-4 | Hartree/Bohr | 最大梯度分量的收敛阈值 |
| `convergence_drms` | f64 | 1.2e-3 | Angstrom | 位移 RMS 的收敛阈值 |
| `convergence_dmax` | f64 | 1.8e-3 | Angstrom | 最大位移分量的收敛阈值 |

### 坐标系统

- `coordsys`：内坐标系统选择，缺省为 `"tric"`。可选值包括 `"cart"`（笛卡尔坐标）、`"prim"`（基元内坐标）、`"dlc"`（去局域内坐标）、`"hdlc"`（混合去局域内坐标）、`"tric"`（Translation-Rotation Internal Coordinates，缺省）以及 `"tric-p"`。对于大多数体系，缺省的 `"tric"` 是最优选择。

### Hessian 与 Trust Radius 参数

以下参数用于控制优化算法的 Hessian 更新策略和步长管理，通常在处理 QM/MM 体系、含噪声梯度的体系或各向异性强的体系时需要调整。

| 关键词 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `hessian` | String | `"never"` | Hessian 矩阵的计算策略（见下方说明） |
| `reset` | bool | `true` (稳态) / `false` (过渡态) | 当 Hessian 特征值低于 `epsilon` 时，是否重置为初始猜测 Hessian |
| `trust` | f64 | 0.1 | 初始 trust radius (Angstrom) |
| `tmax` | f64 | 0.3 | 最大 trust radius (Angstrom) |
| `tmin` | f64 | 1.0e-4 | 最小 trust radius (Angstrom) |
| `epsilon` | f64 | 1.0e-5 | Hessian 重置的特征值阈值，仅当 `reset = true` 时生效 |
| `subfrctor` | i32 | 1 | 投影掉梯度中净力/净力矩分量：0=不投影，1=自动（缺省），2=强制投影 |
| `usedmax` | bool | `false` | 是否用最大位移分量（而非 RMS）判断 trust radius |

`hessian` 关键词的可选值：

| 值 | 说明 |
|---|---|
| `"never"` | 不做 Hessian 计算（缺省，稳态搜索） |
| `"first"` | 仅对初始结构计算 Hessian（缺省，过渡态搜索） |
| `"last"` | 仅对优化后的结构计算 Hessian |
| `"first+last"` | 初始和优化后的结构各计算一次 Hessian |
| `"stop"` | 不做构型优化，仅计算初始结构的 Hessian |
| `"each"` | 每步都计算 Hessian（计算量大，慎用） |

**典型调参场景：**

- **QM/MM 体系**：DFT 梯度中常含微量净力矩噪声，可能导致结构缓慢旋转而力不收敛。设置 `subfrctor = 2` 可强制投影掉净力和净力矩分量。
- **含噪声梯度的体系**：若 BFGS 更新频繁失败（日志中出现 "Eigenvalues below ... returning guess"），可设置 `reset = false` 保留 Hessian 并通过对角 shift 使其正定，而非每次都重置为初始猜测。
- **各向异性体系**：各方向力常数差异较大时，设置 `usedmax = true` 用最大位移分量代替 RMS 判断 trust radius，可获得更稳定的优化行为。

### 过渡态搜索与 IRC

| 关键词 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `transition` | bool | `false` | 设为 `true` 开启过渡态搜索（鞍点优化） |
| `irc` | bool | `false` | 设为 `true` 开启 IRC（内禀反应坐标）路径计算 |
| `irc_direction` | String | `"both"` | IRC 计算方向：`"forward"`、`"backward"` 或 `"both"` |

开启过渡态搜索时，程序会自动将 `hessian` 的缺省值调整为 `"first"`，以在初始结构上计算 Hessian 矩阵，确保优化器能正确识别鞍点方向。

### 频率分析与热化学

| 关键词 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `frequency` | bool | `true` | 当 Hessian 矩阵可用时，是否进行频率计算和热化学分析 |
| `thermo` | [f64; 2] | `[300.0, 1.0]` | 热力学分析的状态参数：[温度 (K), 压强 (bar)] |

### 输出控制

| 关键词 | 类型 | 缺省值 | 说明 |
|---|---|---|---|
| `prefix` | String | `"GeomeTRIC"` | geomeTRIC 生成文件的前缀名 |
| `verbose` | i32 | 0 | 输出详细程度：0=简洁，1=正常，2=详细，3=底层函数输出 |

## L-BFGS 优化器

L-BFGS 引擎是 REST 内置的轻量级优化器，基于 Rust 的 `liblbfgs` 库实现。它始终使用数值梯度，且不支持约束优化（固定原子）和过渡态搜索。对于需要这些功能的场景，请使用缺省的 geomeTRIC 引擎。

相关关键词：

- `nforce_displacement`：数值力计算的有限差分位移步长，单位 Bohr，缺省值为 `0.0013`。

## 完整示例

### 过渡态搜索与频率分析

```toml
[ctrl]
    job_type =                  "opt"
    xc =                        "b3lyp"
    basis_path =                "def2-TZVP"
    auxbas_path =               "def2-universal-jkfit"

[ctrl.geometric_pyo3]
    transition =                true
    hessian =                   "first+last"
    frequency =                 true
    thermo =                    [398.0, 1.5]

[geom]
    name = "ts_example"
    unit = "angstrom"
    position = '''
        ...
    '''
```


