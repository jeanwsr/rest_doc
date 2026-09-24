# TDDFT 模块

REST 的 TDDFT 模块 (`src/ri_tddft/`) 实现了基于 RI 加速的线性响应含时密度泛函理论。该模块复用 `solvers::davidson` 的通用 Davidson 求解器与 `solvers::feast` 围道积分求解器，以及 BSE (`ri_bse`) 模块中的响应方程求解器与偶极工具，但使用 KS 轨道能量（而非 GW 准粒子能量）构建对角项。

模块能力概览：

- **限制性 (RHF/RKS) 参考**：`tddft_spin` 选择单重态 (`"singlet"`)、三重态 (`"triplet"`，仅 AO 模式) 或两者都算 (`"both"`，仅 AO 模式；供 PySOC 导出使用)；
- **非限制 (UHF/UKS) 参考 (UTDDFT)**：`spin_polarization = true`；响应只有单一自旋耦合通道 (α/β 激发通道经自旋无关的库仑核耦合)，`tddft_spin` 不适用；
- **范围分离杂化 (RSH) 泛函**：响应交换按 `c_LR·K_full + (c_SR − c_LR)·K_SR` 分解 (见下文「范围分离杂化」)；
- **SCF 稳定性分析**：`stability`/`check_stab` 关键字，基于 AO 模式 (A+B) 轨道 Hessian (见 [tddft-stability](tddft-stability.md))；
- **激发态解析梯度**：`tddft_grad_state` 关键字，PySCF `grad/tdrks.py` 的 RI 移植 (见 [tddft-grad](tddft-grad.md))。

模块支持两种内核实现模式，由输入关键字 `tddft_mode` 控制：

- **MO 模式** (`tddft_mode = "mo"`，默认)：预先将三中心 RI 积分变换到 MO 基，矩阵-矢量积全部在 MO 振幅空间完成（见 [tddft-mo](tddft-mo.md)）；
- **AO 模式** (`tddft_mode = "ao"`)：Davidson 迭代仍停留在 MO 振幅空间，但每次矩阵-矢量积时构造 AO 跃迁密度，库仑与交换项通过 `ri_jk` 的密度驱动接口计算，XC 核通过 `dft::numint_matmul` 在数值格点上批量计算（见 [tddft-ao](tddft-ao.md)）。

两种模式共享同一数据结构 `TDDFTData` (`src/ri_tddft/tddft.rs`)：模式相关的成员以 `Option` 形式存在；非限制参考下每自旋扇区的成员 (`c_occ`/`c_vir`/`psi_occ`/`ri_terms` 等) 以 `Vec` (每扇区一项) 存储，扇区数由 `TDDFTData::n_sectors()` 给出 (RHF 为 1，UHF 为 2)。交换混合系数不存储在数据中，而是在每次矩阵-矢量积时由 `scf.mol.xc_data` 现场导出。

```{toctree}
---
maxdepth: 1
---

tddft-mo
tddft-ao
tddft-stability
tddft-grad
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
- 对限制性参考的单重态：$\kappa_c = 2$
- 对限制性参考的三重态：$\kappa_c = 0$

非限制参考的库仑项以单位权重 ($\kappa_c = 1$) 耦合 $\alpha$/$\beta$ 两个自旋扇区（见下文「非限制参考」），因此上述按限制性自旋适配的 $\kappa_c$ 因子不再适用。REST 的矩阵-矢量积接口以 `xlet: char` 标记通道：`'S'`/`'T'` 对应限制性单重/三重态，`'R'` 为非自旋适配的通用标记（供非限制参考与稳定性分析的 UHF 路径使用）。

### XC 核的自旋通道

XC 核与自旋通道相关 (CPL, 256, 454)：

$$
f_s = f_{\uparrow\uparrow} + f_{\uparrow\downarrow} \qquad
f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}
$$

其中 $f_{\uparrow\uparrow}$、$f_{\uparrow\downarrow}$ 是自旋分辨的二阶泛函导数。非极化 (unpolarized) 的 libxc 求值只能给出 $f_u = \frac{1}{2}(f_{\uparrow\uparrow} + f_{\uparrow\downarrow})$，因此：

- 单重态核 $f_s = 2 f_u$：直接由非极化求值得到；
- 三重态核 $f_t$：无法由非极化求值获得，必须以自旋极化方式在 $\rho_\uparrow = \rho_\downarrow = \rho/2$ 处求值后，沿反对称方向组合（即 PySCF `nr_rks_fxc_st` 的做法）。REST 目前仅在 **AO 模式**支持三重态；MO 模式下请求三重态会直接报错。
- `tddft_spin = "both"`：依次以单重、三重通道各解一次本征值问题（`tddft_main` 内部的 `run_spin` 闭包按通道分别准备 `TDDFTData`，两通道**不共享**自旋适配核）；输出顺序与 JSON 字段顺序均为先单重后三重。

### 非限制参考 (UTDDFT)

非限制 (UHF/UKS) 参考的响应只有**单一**自旋耦合通道：激发空间为 α、β 两扇区占据→虚轨道旋转的拼接 $[z_\alpha; z_\beta]$，两扇区通过自旋无关的库仑核耦合（$J[\sum_\tau z^\tau]$），不存在限制性形式中「因子 2 / 因子 0」的自旋适配对，因此 `tddft_spin` 不适用（显式给出即报错）。XC 核为自旋分辨核 $f_{\sigma_1\sigma_2}[g,\alpha,\beta]$（无单重/三重因子），MO 模式存于 `fxc_u` (spin-resolved `FXCMatvecDataUnrestricted`)，AO 模式存于自旋极化的 `fxc_eff: [n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2]`。虚轨道截断 (`tddft_cutoff_energy`) 在 α/β 通道独立解析（`tddft_occupation_parameters_u`）；冻结芯（`mol.start_mo`，由 `frozen_core_postscf` 控制）为两通道共享；空扇区 (如 β 无占据) 以零维扇区参与。振幅后处理 (跃迁偶极、振子强度、主导跃迁打印) 遵循 PySCF `uhf.py` 约定。

MO 模式与 AO 模式均支持非限制参考。求解器分层与限制性情形相同，稠密对角化阈值独立放宽：MO-U TDA $\dim \le 15$ (`dsyev`)、MO-U Full LR $\dim \le 80$ (直接构造非厄米 $[\mathbf{A}\ \mathbf{B};-\mathbf{B}\ -\mathbf{A}]$ 并以 `dgeev` 对角化)。

### 范围分离杂化 (RSH) 泛函

对范围分离杂化泛函，HF 交换算符按 $1/r_{12} = \mathrm{erfc}(\omega r_{12})/r_{12} + \mathrm{erf}(\omega r_{12})/r_{12}$ 分解为短程 (SR) 与长程 (LR) 两部分。REST 的响应交换约定 (与 `scf_io` 基态 Fock 构建一致) 为

$$
K^{\text{resp}} = c_{LR}\, K_{\text{full}} + (c_{SR} - c_{LR})\, K_{\text{SR}}
$$

其中 $K_{\text{full}}$ 为普通 $1/r_{12}$ 交换、$K_{\text{SR}}$ 为 $\mathrm{erfc}(\omega r_{12})/r_{12}$ 交换，$c_{LR}$、$c_{SR}$ 为泛函的 RSH 参数 (`xc_data.rsh_params()` → $(\omega, c_{LR}, c_{SR})$)。该系数**不预存**在 `TDDFTData` 中，而是每次矩阵-矢量积时由 `scf.mol.xc_data` 现场导出；仅当 $|c_{SR} - c_{LR}| > 10^{-12}$ 时才构建 SR 三中心积分张量 (MO 模式为 `RITensorTerms` 的 `oo_sr`/`vv_sr`/`ov_sr` 三元组；AO 模式要求 `scf.rimatr_sr` 存在，即由同泛函的 SCF 建立)。对普通杂化泛函，$c_{SR} - c_{LR} = 0$，公式退化为 $c_x K_{\text{full}}$。

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
| `tddft_method` | String | `"lr"` | `"tda"` (Tamm-Dancoff 近似) 或 `"lr"` (全线性响应) |
| `tddft_mode` | String | `"mo"` | `"mo"` (MO 基 RI 张量) 或 `"ao"` (AO 跃迁密度核) |
| `tddft_spin` | String | `"singlet"` | `"singlet"` / `"triplet"` / `"both"`；三重态与 both 仅 AO 模式、且仅限制性参考支持（非限制参考显式给出即报错） |
| `nroots` | Integer | 6 | 求解的激发态数目 |
| `tddft_use_optimized_fxc` | Bool | `true` | 是否使用 rayon 并行的 fxc 矩阵-矢量积核 |
| `davidson_tol` | Float | `1e-10` | Davidson 收敛阈值：$\|r\| < \sqrt{\varepsilon}$ 且 $|\Delta E| < \varepsilon$ |
| `davidson_max_iter` | Integer | 50 | Davidson 最大迭代次数 |
| `davidson_max_subspace` | Integer | 60 | 子空间容量上限（实际子空间维度 = `max(4 × nroots, davidson_max_subspace)`，并以激发空间维度截断；Full LR 需要较大子空间在首次重启前收敛，不建议调小） |
| `grid_batch` | Bool | `true` | 仅 AO 模式：XC 核求值按格点分批，避免完整 AO-on-grid 张量常驻内存；MO 模式下忽略 |
| `tddft_ao_rik_driver` | String | `"semitrans"` | 仅 AO 模式：交换 K 驱动方式——`"semitrans"` (占据侧半转换收缩，默认)、`"dm"` (精确批量)、`"lowrank"` (逐向量 SVD 低秩) |
| `tddft_fxc_driver` | String | `"semitrans"` | 仅 AO 模式：fxc 驱动方式——`"semitrans"` (C_vir 折入振幅，虚轨道侧直接与格点裸 AO 收缩，无需形成 psi_vir) 或 `"mo"` (缓存占据侧格点投影 + 虚轨道侧流式，MO 模式 fxc 算法) 或 `"dm"` (组装密度 NIMatmul 回退路径；未知取值告警并回退到此) |
| `tddft_svd_tol` | Float | `1e-6` | 仅 AO 模式：低秩 K 的相对奇异值阈值（保留 $\sigma_i \ge \varepsilon \sigma_{\max}$） |
| `stability` | String | `"off"` | SCF 稳定性分析 (`"internal"` / `"external"` / `"full"` / `"auto"`)；非 `"off"` 时本任务只做稳定性分析（见 [tddft-stability](tddft-stability.md)） |
| `stability_nroots` | Integer | 3 | 稳定性 Hessian 求解的最低本征值数目 |
| `stability_tol` | Float | `1e-8` | 稳定性 Hessian 的 Davidson 收敛阈值 |
| `tddft_grad_state` | Integer | 0 | 激发态解析梯度 (1-based，0 关闭；仅限制性参考，见 [tddft-grad](tddft-grad.md)) |
| `pysoc` | Bool | `false` | 导出 PySOC JSON (`rest_pysoc_export.json`)；要求限制性参考 + `tddft_spin = "both"` |
| `tddft_cutoff_energy` | Float | `1e6` | 虚轨道能量截断 (Hartree)；非限制下 α/β 通道独立解析 |
| `tddft_feast_solver` | Bool | `false` | FEAST 围道积分求解器（仅 MO 模式、限制性参考） |
| `response_tddft` | Bool | `false` | 启用响应 TDDFT（MO 实现、限制性参考） |
| `response_tddft_solver` | String | `"klopper"` | `"pople"` / `"gmres"` / `"klopper"` / `"dense"` |

NOTE: 顶层 `[ctrl]` 区块的 `check_stab` 关键字与 `[tddft] stability` 接受相同取值；两者同时给出时后者优先。

## 计算流程

### 本征值 TDDFT (`tddft_main`)

```
tddft_main(scf)
    │
    ├── Step 1: 读取 TDDFT 控制参数 + 参考类型门控
    │   ├── tddft_method (tda/lr), tddft_spin, tddft_mode, nroots, ...
    │   ├── ROHF 参考 → 报错 (请改用 spin_polarization = true)
    │   ├── 非限制参考 → tddft_spin 必须缺省；pysoc 需限制性 + "both"
    │   └── 限制性参考 → 解析 tddft_spin (singlet/triplet/both) → xlet
    │       三重态/both 且 MO 模式 → 报错 (AO 模式限定)
    │
    ├── Step 2: 确定轨道扇区
    │   ├── 限制性: tddft_occupation_parameters() → (start_mo, occ_size, vir_size, dim)
    │   └── 非限制: tddft_occupation_parameters_u() → [α 扇区, β 扇区]
    │       虚轨道截断 (tddft_cutoff_energy) 在此处理；冻结芯为 mol.start_mo (由 frozen_core_postscf 控制)
    │
    ├── Step 3: 生成初始猜测 + 对角预条件器
    │   └── build_hdiag() + generate_initial_guess() (来自 solvers/davidson)
    │
    ├── Step 4: 按自旋通道求解 (run_spin 闭包，每通道独立准备数据)
    │   ├── 数据准备 (每通道独立，自旋适配核不共享):
    │   │   ├── MO 模式: prepare_mo_data() → fxc/fxc_u 表 + 每扇区 MO 基 RI 束 (ri_terms)
    │   │   └── AO 模式: prepare_ao_data_with_spin(scf, Some(spin)) → c_occ/c_vir + NIMatmul + 原始 XC 核
    │   │
    │   └── 求解器分层分派:
    │       ├── FEAST (仅限制性 + MO): feast_solve_tddft_tda/lr()
    │       ├── 稠密对角化 (小系统):
    │       │   ├── 限制性 dim ≤ 15: build_a/build_b → dsyev (LR 经 (A−B) 对称化约化，失败回退 TDA)
    │       │   ├── MO-U TDA dim ≤ 15: dsyev
    │       │   └── MO-U LR dim ≤ 80: 显式 [A B; -B -A] + dgeev
    │       └── Davidson 迭代:
    │           ├── AO 模式: tda/lr_davidson_solver_batched (批量闭包)
    │           ├── MO 限制性: tda/lr_davidson_solver (逐向量闭包)
    │           └── MO 非限制: 逐向量 Davidson (拼接 [α;β] 振幅)
    │
    ├── Step 5: 逐通道输出激发能与振子强度
    │   ├── transition_dipole_square(), normalize() (非限制用 *_u 变体，PySCF uhf.py 约定)
    │   └── tddft_spin = "both": 先单重后三重；pysoc = true 时导出 rest_pysoc_export.json
    │
    └── 汇总写入 rest_results.json 的 "tddft" 字段 (energies, oscillator_strength)
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

### SCF 稳定性分析 (`stability::stability`)

```
main_driver (SCF 收敛后、激发态计算之前)
    │
    ├── 解析稳定性模式: [tddft] stability 优先，[ctrl] check_stab 兜底
    │
    ├── stability::stability(scf, mode)
    │   ├── internal: RHF/RKS → H = 4(A^S+B^S)；UHF/UKS → H = 2(A+B) (拼接 [α;β])
    │   ├── external (仅 RHF/RKS): H = (A^T+B^T) (RHF→UHF 方向)
    │   └── batched Davidson 求最低 nroots 个本征值 → λ_min < -1e-5 即不稳定
    │
    └── 结果打印 + rest_results.json 的 "stability" 字段 (internal/external roots)
```

详见 [tddft-stability](tddft-stability.md)。

### 激发态解析梯度 (`tddft_grad::TddftGradEngine`)

`tddft_grad_state > 0` 时，梯度任务在 `main_driver::eval_force` 中把激发态响应梯度叠加到基态梯度上。模块细节待补充，见 [tddft-grad](tddft-grad.md)。

## 代码结构

`src/ri_tddft/` 目录包含：

| 文件 | 职责 |
|------|------|
| `mod.rs` | 子模块声明与公开 API 导出 (`tddft_main`, `response_tddft`, `TddftGradEngine`) |
| `tddft.rs` | 共享数据结构 `TDDFTData`、`TDDFTMode`/`FxcDriver` 枚举，与两个构造器 (`prepare_mo_data`, `prepare_ao_data_with_spin`)，以及稠密小系统路径的模式分派构造器 `build_a`/`build_b` |
| `tddft_solver.rs` | 本征值 TDDFT 总调度器 (`tddft_main`)。协调参考类型门控、扇区解析、逐自旋通道数据准备与求解器分派 (FEAST/稠密/Davidson)、结果打印与 JSON 汇总 |
| `matvec.rs` | MO 模式矩阵-矢量积。实现 A 块和 B 块的矩阵-矢量积 (`a_matvec`, `b_matvec`，含 RSH 双系数交换)，每扇区 RI 束 `RITensorTerms` |
| `matvec_ao.rs` | AO 模式矩阵-矢量积。跃迁密度构造、`ri_jk` 批量 J/K、`numint_matmul` 批量 fxc、批量与稠密两条路径 |
| `response.rs` | 响应 TDDFT 求解器 (`response_tddft`)。实现 Pople、GMRES、Klopper、稠密 LU 四种求解后端，以及极化密度导出 |
| `stability.rs` | SCF 波函数稳定性分析 (`stability`, `StabilityReport`)。AO 模式 (A+B) 轨道 Hessian + batched Davidson |
| `tddft_grad.rs` | 激发态解析梯度 (`TddftGradEngine`)。PySCF `grad/tdrks.py`/`tdrhf.py` 的 RI 移植 |
| `utils.rs` | 工具函数。`tddft_occupation_parameters`/`_u`/`tddft_sector_params` (轨道扇区计算)、`tddft_get_submatrix` (RI 子矩阵抽取)、`compute_tddft_dipole_matrix` (跃迁偶极) |
| `feast_solver.rs` | FEAST 求解器封装 (仅 MO 模式、限制性参考)。将 TDDFT 的矩阵-矢量积适配到通用 `solvers::feast` 接口 |

输入参数定义位于 `src/ctrl_io/tddft_parameters.rs`。

## 集成点

### 上游依赖

| 依赖模块 | 用途 |
|----------|------|
| `scf_io::SCF` | 提供 MO 系数、KS 轨道能量、RI 积分 (`rimatr`，RSH 另有 `rimatr_sr`)、数值格点、分子信息等核心数据 |
| `dft::num_int` | 提供 `FXCMatvecData`/`FXCMatvecDataUnrestricted` (MO 模式 XC 核数据) 和 `fxc_matvec()` (XC 核矩阵-矢量积) |
| `dft::numint_matmul` | AO 模式 XC 核：`NIMatmul` (格点 AO 缓存、批量密度构造与核收缩)、`eval_vxc_fxc_from_rho` (原始 XC 核) |
| `dft::xceff` | libxc 求值封装 (`libxc_eval_eff`, `determine_den_type`)；三重态自旋极化核在此之上组合 |
| `ri_jk` | AO 模式 J/K：`get_vj_ri_incore_nonsym`、`get_vk_ri_incore_dm`、`get_vk_ri_incore_dm_lowrank`、`get_vk_ri_incore_coeff_pair` (记号见 [ri-jk 文档](../ri-jk/index.md)) |
| `ri_bse` | 提供 Coulomb 贡献 (`coulomb_contribution`)、响应方程求解器 (Pople/GMRES/Klopper，`ri_bse::response` 的适配复用)、偶极工具 (`dipoles::normalize` 等)、`pysoc_export` |
| `solvers::davidson` | 提供通用 Davidson 求解器：逐向量接口 (`davidson_solver`, `lr_davidson_solver`) 与批量接口 (`davidson_solver_batched`, `lr_davidson_solver_batched`) |
| `solvers::feast` | 提供通用 FEAST 围道积分算法 |
| `ri_cphf` | 梯度模块的 CPHF Z-vector 求解 (`CPHFSolverPySCF`) |
| `ri_gw::gw_grad` | 梯度模块复用的 RI 原始张量与原子导数分块 (`RawRiTensors`, `AtomDerivBlocks`) |

### 与 BSE 的关系

TDDFT 模块重用了 `ri_bse` 模块的响应求解器与偶极工具等基础设施，但有两项关键差异：

1. **对角项**：TDDFT 使用 KS 轨道能量差 $\varepsilon_a - \varepsilon_i$，而 BSE 使用 GW 准粒子能量差。
2. **数据流**：TDDFT 位于 `main_driver` 的直接调用链上（与 PT2、RPA 并列），而 BSE 需要先运行 GW 计算。
