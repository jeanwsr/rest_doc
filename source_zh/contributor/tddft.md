# TDDFT 模块

REST 的 TDDFT 模块 (`src/ri_tddft/`) 实现了基于 RI 加速的线性响应含时密度泛函理论。该模块复用了 BSE (`ri_bse`) 模块中的 Davidson 求解器、响应方程求解器和 FEAST 围道积分求解器等基础设施，但使用 KS 轨道能量（而非 GW 准粒子能量）构建对角项。

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
A_{ia,jb} = (\varepsilon_a - \varepsilon_i) \delta_{ij} \delta_{ab} + 2(ia|jb) - c_x (ij|ab) + f_{ia,jb}^{\text{xc}}
$$

$$
B_{ia,jb} = 2(ia|bj) - c_x (ib|aj) + f_{ia,jb}^{\text{xc}}
$$

这里 $\varepsilon_i, \varepsilon_a$ 为 KS 轨道能量（占据和虚轨道），$(ia|jb)$ 为 MO 基下的双电子积分，$f^{\text{xc}}$ 为 XC 核 (exchange-correlation kernel)，$c_x$ 为 HF 交换混合系数。

- 对纯泛函：$c_x = 0$
- 对杂化泛函：$c_x = \alpha_{\text{hybrid}}$
- 对单重态：库仑因子为 2
- 对三重态：库仑因子为 0

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

### RI 加速的矩阵-矢量积

所有矩阵-矢量积均通过 RI 分解加速。以库仑项为例：

$$
\sum_{jb} (ia|jb) z_{jb} = \sum_{P} B_{ia}^P \sum_{jb} B_{jb}^P z_{jb}
$$

其中 $B_{ia}^P$ 为三中心 RI 积分。交换项通过 DGEMM 操作实现：

$$
K_A \cdot z = -\alpha \cdot \mathbf{RI}_{OO}^T \cdot (\mathbf{RI}_{VV} \cdot z^T)^T
$$

## 计算流程

### 本征值 TDDFT (`tddft_main`)

```
tddft_main(scf)
    │
    ├── Step 1: 读取 TDDFT 控制参数
    │   └── tddft_method (tda/lr), tddft_spin (singlet/triplet), nroots, ...
    │
    ├── Step 2: 确定轨道维度
    │   └── tddft_occupation_parameters() → (start_mo, occ_size, vir_size, dim)
    │       冻结芯 (< -2.0 Ha) 和虚轨道截断 (tddft_cutoff_energy) 在此处理
    │
    ├── Step 3: 准备 fxc 数据
    │   └── prepare_fxc_data() → FXCMatvecData (含 alpha_hybrid)
    │
    ├── Step 4: 获取并重塑 RI 积分
    │   └── tddft_get_submatrix() → ri_ov, ri_oo, ri_vv
    │       reshape → ri_oo_exch, ri_vv_exch, ri_ov_exch (交换项用)
    │
    ├── Step 5: 构建对角预条件器
    │   └── matvec::build_hdiag() → hdiag = ε_a - ε_i
    │
    ├── Step 6: 生成初始猜测
    │   └── generate_initial_guess() (来自 solvers/davidson)
    │
    ├── Step 7: 分派求解器
    │   ├── dim ≤ 15: 稠密对角化 (LAPACK dsyev/dsyevd)
    │   ├── FEAST 启用: feast_solve_tddft_tda() 或 feast_solve_tddft_lr()
    │   └── 默认: Davidson 迭代求解器
    │       ├── TDA: tda_davidson_solver()
    │       └── Full LR: lr_davidson_solver()
    │
    └── Step 8: 输出激发能与振子强度
        └── transition_dipole_square(), normalize()
```

### 响应 TDDFT (`response_tddft`)

```
response_tddft(scf)
    │
    ├── 准备 fxc 数据 + RI 积分 (同上)
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

`src/ri_tddft/` 目录包含 5 个源文件和 1 个模块声明文件：

| 文件 | 职责 |
|------|------|
| `mod.rs` | 子模块声明与公开 API 导出 (`tddft_main`, `response_tddft`) |
| `tddft_solver.rs` | 本征值 TDDFT 总调度器 (`tddft_main`)。协调参数解析、积分准备、求解器分派和输出打印 |
| `response.rs` | 响应 TDDFT 求解器 (`response_tddft`)。实现 Pople、GMRES、Klopper、稠密 LU 四种求解后端，以及极化密度导出 |
| `matvec.rs` | 矩阵-矢量积核心。实现 A 块和 B 块的矩阵-矢量积 (`a_matvec`, `b_matvec`)，包含库仑、精确交换和 fxc 三部分贡献 |
| `utils.rs` | 工具函数。`tddft_occupation_parameters` (轨道维度计算)、`tddft_get_submatrix` (RI 子矩阵抽取)、`compute_tddft_dipole_matrix` (跃迁偶极) |
| `feast_solver.rs` | FEAST 求解器封装。将 TDDFT 的矩阵-矢量积适配到通用 `ri_bse::feast_solver::feast()` 接口 |

输入参数定义位于 `src/ctrl_io/tddft_parameters.rs`。

## 集成点

### 上游依赖

| 依赖模块 | 用途 |
|----------|------|
| `scf_io::SCF` | 提供 MO 系数、KS 轨道能量、RI 积分、分子信息等核心数据 |
| `dft::num_int` | 提供 `FXCMatvecData` (XC 核数据) 和 `fxc_matvec()` (XC 核矩阵-矢量积) |
| `ri_bse` | 提供 Coulomb 贡献 (`coulomb_contribution`)、响应方程求解器 (Pople/GMRES/Klopper)、偶极积分工具、FEAST 算法 |
| `solvers::davidson` | 提供通用 Davidson 迭代求解器 (`tda_davidson_solver`, `lr_davidson_solver`) |

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
