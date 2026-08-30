# TDDFT 模块

REST 的 TDDFT 模块 (`src/ri_tddft/`) 实现了基于 RI 加速的线性响应含时密度泛函理论。该模块复用了 BSE (`ri_bse`) 模块中的 Davidson 求解器、响应方程求解器和 FEAST 围道积分求解器等基础设施，但使用 KS 轨道能量（而非 GW 准粒子能量）构建对角项。

模块支持两种内核实现模式，由输入关键字 `tddft_mode` 控制：

- **MO 模式** (`tddft_mode = "mo"`，默认)：预先将三中心 RI 积分变换到 MO 基，矩阵-矢量积全部在 MO 振幅空间完成（见 [tddft-mo](tddft-mo.md)）；
- **AO 模式** (`tddft_mode = "ao"`)：Davidson 迭代仍停留在 MO 振幅空间，但每次矩阵-矢量积时构造 AO 过渡密度，库仑与交换项通过 `ri_jk` 的密度驱动接口计算，XC 核通过 `dft::numint_matmul` 在数值格点上批量计算（见 [tddft-ao](tddft-ao.md)）。

两种模式共享同一数据结构 `TDDFTData` (`src/ri_tddft/tddft.rs`)：模式相关的成员以 `Option` 形式存在，公共标量（`alpha_hybrid`）直接存储。

```{toctree}
---
maxdepth: 1
---

tddft-mo
tddft-ao
```

## 理论背景

### Casida 方程

在 KS-DFT 框架下，线性响应 TDDFT 的本征值问题可写为如下非厄米 Casida 方程：

$$
\begin{pmatrix}
\mathbf{A} & \mathbf{B} \\
-\mathbf{B} & -\mathbf{A}
\end{pmatrix}
\begin{pmatrix}
\mathbf{X} \\
\mathbf{Y}
\end{pmatrix}
= \omega
\begin{pmatrix}
\mathbf{X} \\
\mathbf{Y}
\end{pmatrix}
$$

其中矩阵块 $\mathbf{A}$ 和 $\mathbf{B}$ 的矩阵元 (在 MO 基下) 为：

$$
A_{ia,jb} = (\varepsilon_a - \varepsilon_i) \delta_{ij} \delta_{ab} + \kappa_c (ia|jb) - c_x (ij|ab) + f_{ia,jb}^{\text{xc}}
$$

$$
B_{ia,jb} = \kappa_c (ia|bj) - c_x (ib|aj) + f_{ia,jb}^{\text{xc}}
$$

这里 $\varepsilon_i, \varepsilon_a$ 为 KS 轨道能量（占据和虚轨道），$(ia|jb)$ 为 MO 基下的双电子积分，$f^{\text{xc}}$ 为 XC 核 (exchange-correlation kernel)，$c_x$ 为 HF 交换混合系数，$\kappa_c$ 为库仑耦合因子：

- 对纯泛函：$c_x = 0$
- 对杂化泛函：$c_x = \alpha_{\text{hybrid}}$
- 对单重态：$\kappa_c = 2$
- 对非极化 ('R') 情形：$\kappa_c = 1$
- 对三重态：$\kappa_c = 0$

### XC 核的自旋通道

XC 核与自旋通道相关 (CPL, 256, 454)：

$$
f_s = f_{\uparrow\uparrow} + f_{\uparrow\downarrow} \qquad
f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}
$$

其中 $f_{\uparrow\uparrow}$、$f_{\uparrow\downarrow}$ 是自旋分辨的二阶泛函导数。非极化 (unpolarized) 的 libxc 求值只能给出 $f_u = \frac{1}{2}(f_{\uparrow\uparrow} + f_{\uparrow\downarrow})$，因此：

- 单重态核 $f_s = 2 f_u$：直接由非极化求值得到；
- 三重态核 $f_t$：无法由非极化求值获得，必须以自旋极化方式在 $\rho_\uparrow = \rho_\downarrow = \rho/2$ 处求值后，沿反对称方向组合（即 PySCF `nr_rks_fxc_st` 的做法）。REST 目前仅在 **AO 模式**支持三重态；MO 模式下请求三重态会直接报错。

### Tamm-Dancoff 近似 (TDA)

忽略 $\mathbf{B}$ 矩阵块 (即令 $\mathbf{Y} = \mathbf{0}$)，简化方程至标准厄米本征值问题：

$$
\mathbf{A} \mathbf{X} = \omega \mathbf{X}
$$

TDA 通常给出与实验吻合较好的激发能，且求解效率更高（无需处理非厄米矩阵）。

### 全线性响应 (Full LR)

对于全 Casida 方程，可通过如下变换转化为厄米形式：

$$
(\mathbf{A} - \mathbf{B})^{1/2}(\mathbf{A} + \mathbf{B})(\mathbf{A} - \mathbf{B})^{1/2} \mathbf{Z} = \omega^2 \mathbf{Z}
$$

求解后通过反变换 $\mathbf{X} = (\mathbf{A} - \mathbf{B})^{1/2} \mathbf{Z} / \sqrt{\omega}$ 和 $\mathbf{Y} = \mathbf{X} - \mathbf{Z} / \sqrt{\omega}$ 恢复 X, Y 矢量。

### 响应 TDDFT

频率域响应 TDDFT 求解如下 4 分量非厄米线性系统：

$$
\begin{pmatrix}
\mathbf{A} - \omega & \mathbf{B} & -\gamma & \mathbf{0} \\
\mathbf{B} & \mathbf{A} + \omega & \mathbf{0} & \gamma \\
\gamma & \mathbf{0} & \mathbf{A} - \omega & \mathbf{B} \\
\mathbf{0} & -\gamma & \mathbf{B} & \mathbf{A} + \omega
\end{pmatrix}
\begin{pmatrix}
\mathbf{r}^+ \\
\mathbf{r}^- \\
\mathbf{i}^+ \\
\mathbf{i}^-
\end{pmatrix}
=
\begin{pmatrix}
\boldsymbol{\mu}_z \\
\boldsymbol{\mu}_z \\
\mathbf{0} \\
\mathbf{0}
\end{pmatrix}
$$

其中 $\omega$ 为外场频率，$\gamma$ 为寿命展宽，$\boldsymbol{\mu}_z$ 为 z 方向跃迁偶极矩矢量。求解后计算频域极化率 $\alpha_{zz}(\omega)$。

响应 TDDFT (`response_tddft`) 目前只有 MO 实现，且与 `tddft_mode` 无关——无论 `tddft_mode` 取何值，响应求解器总是使用 MO 机制。

## 输入关键字

与 TDDFT 相关的主要输入关键字（定义于 `src/ctrl_io/tddft_parameters.rs`）：

| 关键字 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `tddft_method` | String | `"tda"` | `"tda"` (Tamm-Dancoff 近似) 或 `"lr"` (全线性响应) |
| `tddft_mode` | String | `"mo"` | `"mo"` (MO 基 RI 张量) 或 `"ao"` (AO 过渡密度核) |
| `tddft_spin` | String | `"singlet"` | `"singlet"` 或 `"triplet"`；三重态目前仅 AO 模式支持 |
| `nroots` | Integer | 1 | 求解的激发态数目 |
| `grid_batch` | Bool | `true` | 仅 AO 模式：XC 核求值按格点分批，避免完整 AO-on-grid 张量常驻内存；MO 模式下忽略 |
| `tddft_ao_rik_driver` | String | `"semitrans"` | 仅 AO 模式：交换 K 驱动方式——`"semitrans"` (占据侧半变换收缩，默认)、`"dm"` (精确批量)、`"lowrank"` (逐向量 SVD 低秩) |
| `tddft_fxc_driver` | String | `"mo"` | 仅 AO 模式：fxc 驱动方式——`"mo"` (缓存占据侧格点投影 + 虚轨道侧流式，MO 模式 fxc 算法，occ/vir 约化收缩) 或 `"dm"` (组装密度 NIMatmul 回退路径) |
| `tddft_svd_tol` | Float | `1e-6` | 仅 AO 模式：低秩 K 的相对奇异值阈值（保留 $\sigma_i \ge \varepsilon \sigma_{\max}$） |
| `tddft_feast_solver` | Bool | `false` | FEAST 围道积分求解器（仅 MO 模式） |
| `response_tddft` | Bool | `false` | 启用响应 TDDFT（MO 实现） |
| `response_tddft_solver` | String | `"klopper"` | `"pople"` / `"gmres"` / `"klopper"` / `"dense"` |

## 计算流程

### 本征值 TDDFT (`tddft_main`)

```
tddft_main(scf)
    │
    ├── Step 1: 读取 TDDFT 控制参数
    │   └── tddft_method (tda/lr), tddft_spin, tddft_mode, nroots, ...
    │
    ├── Step 2: 确定轨道维度
    │   └── tddft_occupation_parameters() → (start_mo, occ_size, vir_size, dim)
    │       冻结芯 (< -2.0 Ha) 和虚轨道截断 (tddft_cutoff_energy) 在此处理
    │
    ├── Step 3: 准备共享数据 TDDFTData
    │   ├── MO 模式: prepare_mo_data() → fxc 表 + 四个 MO 基 RI 张量
    │   └── AO 模式: prepare_ao_data() → c_occ/c_vir + NIMatmul + 原始 fxc 核表
    │
    ├── Step 4: 构建对角预条件器
    │   └── matvec::build_hdiag() → hdiag = ε_a - ε_i
    │
    ├── Step 5: 生成初始猜测
    │   └── generate_initial_guess() (来自 solvers/davidson)
    │
    ├── Step 6: 分派求解器
    │   ├── dim ≤ 15: 稠密对角化 (build_a/build_b + LAPACK dsyev)
    │   ├── FEAST 启用 (仅 MO): feast_solve_tddft_tda() 或 feast_solve_tddft_lr()
    │   └── 默认: Davidson 迭代求解器 (批量接口)
    │       ├── TDA: davidson_solver_batched()
    │       └── Full LR: lr_davidson_solver_batched()
    │
    └── Step 7: 输出激发能与振子强度
        └── transition_dipole_square(), normalize()
```

### 响应 TDDFT (`response_tddft`)

```
response_tddft(scf)
    │
    ├── 准备 MO 模式数据 (prepare_mo_data，与 tddft_mode 无关)
    │
    ├── 构建 KS 能量对角矢量 + 偶极矢量
    │
    ├── 构建 4 分量矩阵-矢量积闭包
    │   (A, B, fxc, exchange 全部打包为 H_4c * x)
    │
    ├── 选择求解器
    │   ├── "pople":   Pople-Krylov 数值技巧
    │   ├── "gmres":   4 分量 GMRES (重启=30)
    │   ├── "klopper": Klopper 子空间求解器 (默认)
    │   └── "dense":   稠密 LU 分解 (LAPACK)
    │
    ├── 计算动态极化率 α_zz(ω)
    │
    └── (可选) 导出极化密度到空间格点文件
```

## 代码结构

`src/ri_tddft/` 目录包含：

| 文件 | 职责 |
|------|------|
| `mod.rs` | 子模块声明与公开 API 导出 (`tddft_main`, `response_tddft`) |
| `tddft.rs` | 共享数据结构 `TDDFTData` 与两个构造器 (`prepare_mo_data`, `prepare_ao_data`)，以及稠密小系统路径的模式分派构造器 `build_a`/`build_b` |
| `tddft_solver.rs` | 本征值 TDDFT 总调度器 (`tddft_main`)。协调参数解析、数据准备、求解器分派和输出打印 |
| `matvec.rs` | MO 模式矩阵-矢量积。实现 A 块和 B 块的矩阵-矢量积 (`a_matvec`, `b_matvec`) |
| `matvec_ao.rs` | AO 模式矩阵-矢量积。过渡密度构造、`ri_jk` 批量 J/K、`numint_matmul` 批量 fxc、批量与稠密两条路径 |
| `response.rs` | 响应 TDDFT 求解器 (`response_tddft`)。实现 Pople、GMRES、Klopper、稠密 LU 四种求解后端，以及极化密度导出 |
| `utils.rs` | 工具函数。`tddft_occupation_parameters` (轨道维度计算)、`tddft_get_submatrix` (RI 子矩阵抽取)、`compute_tddft_dipole_matrix` (跃迁偶极) |
| `feast_solver.rs` | FEAST 求解器封装 (仅 MO 模式)。将 TDDFT 的矩阵-矢量积适配到通用 `ri_bse::feast_solver::feast()` 接口 |

输入参数定义位于 `src/ctrl_io/tddft_parameters.rs`。

## 集成点

### 上游依赖

| 依赖模块 | 用途 |
|----------|------|
| `scf_io::SCF` | 提供 MO 系数、KS 轨道能量、RI 积分 (`rimatr`)、数值格点、分子信息等核心数据 |
| `dft::num_int` | 提供 `FXCMatvecData` (MO 模式 XC 核数据) 和 `fxc_matvec()` (XC 核矩阵-矢量积) |
| `dft::numint_matmul` | AO 模式 XC 核：`NIMatmul` (格点 AO 缓存、批量密度构造与核收缩) |
| `ri_jk` | AO 模式 J/K：`get_vj_ri_incore_nonsym`、`get_vk_ri_incore_dm`、`get_vk_ri_incore_dm_lowrank` (记号见 [ri-jk 文档](../ri-jk/index.md)) |
| `ri_bse` | 提供 Coulomb 贡献 (`coulomb_contribution`)、响应方程求解器 (Pople/GMRES/Klopper)、偶极积分工具、FEAST 算法 |
| `solvers::davidson` | 提供通用 Davidson 求解器：逐向量接口 (`davidson_solver`, `lr_davidson_solver`) 与批量接口 (`davidson_solver_batched`, `lr_davidson_solver_batched`) |

### 下游消费者

| 调用位置 | 用途 |
|----------|------|
| `main_driver.rs` | 在 SCF 收敛后调用 `tddft_main()` 或 `response_tddft()` |
| `ri_cphf/cphf_solver_pyscf.rs` | 使用 `tddft_occupation_parameters()` 获取轨道维度 |
| `dft/num_int.rs` | 使用 `tddft_occupation_parameters()` 获取轨道维度 |
| `dft/response.rs` | 使用 `tddft_occupation_parameters()` 获取轨道维度 |

### 与 BSE 的关系

TDDFT 模块重用了 `ri_bse` 模块的大量基础设施（求解器、偶极工具、响应求解器），但有两项关键差异：

1. **对角项**：TDDFT 使用 KS 轨道能量差 $\varepsilon_a - \varepsilon_i$，而 BSE 使用 GW 准粒子能量差。
2. **数据流**：TDDFT 位于 `main_driver` 的直接调用链上（与 PT2、RPA 并列），而 BSE 需要先运行 GW 计算。
