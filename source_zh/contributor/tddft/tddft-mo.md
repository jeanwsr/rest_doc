# TDDFT MO 模式：矩阵-矢量积

本页描述 `tddft_mode = "mo"` 下的矩阵-矢量积实现。该模式沿用 RI-JK 的记号体系（见 [RI-JK Incore 算法](../ri-jk/ri-jk-incore.md)），核心思想是**预先将三中心 RI 积分变换到 MO 基**，此后所有矩阵-矢量积都只是 MO 基张量上的 DGEMM/DGEMV 缩并。

## 记号约定

除 [ri-jk 文档](../ri-jk/index.md)的约定外，本页额外使用：

- 占据 MO 指标 $i, j$，虚轨道 MO 指标 $a, b$，辅助基指标 $P, Q$；
- 激发振幅 $z_{ia}^{\mathbb{A}}$，集合指标 $\mathbb{A}$ 遍历 Davidson 子空间中的试探向量（逐向量接口中 $n_\mathrm{set} = 1$）；
- MO 基三中心积分 $B_{ia, P} = \sum_{\mu\nu} C_{\mu i} Y_{\mu\nu, P} C_{\nu a}$，度量已吸收进 $Y$；
- 维度大小 $n_\mathrm{occ}$、$n_\mathrm{vir}$、$n_\mathrm{aux}$、$n_\mathrm{basis}$、$n_\mathrm{grid}$、$n_\mathrm{var}$（密度变量数，LDA 为 1，GGA 为 4）。

下划线角标表示计算时被分批的指标。

## 函数 `prepare_mo_data`

函数路径：`ri_tddft::tddft::prepare_mo_data`

构造 `TDDFTData`：调用 `prepare_fxc_data` 得到 fxc 核表，并经 `tddft_get_submatrix` 从 `scf.rimatr` 抽取三个 MO 基 RI 子张量、按交换项的缩并需求重塑出另外三个：

$$
\begin{aligned}
B_{ia, P} &\;\leftarrow\; \texttt{tddft\_get\_submatrix}(\texttt{'O'}, \texttt{'V'}) && \text{(eq.1)}
\end{aligned}
$$

| 变量名 (`TDDFTData` 成员) | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `ri_ov` | $B_{ia, P}$ | $(P, ia)$ | $(n_\mathrm{aux}, n_\mathrm{occ} n_\mathrm{vir})$ | Coulomb 用 |
| `ri_oo_exch` | $B_{ij, P}$ | $(jP, i)$ | $(n_\mathrm{occ} n_\mathrm{aux}, n_\mathrm{occ})$ | A 块交换用 |
| `ri_vv_exch` | $B_{ab, P}$ | $(P a, b)$ | $(n_\mathrm{aux} n_\mathrm{vir}, n_\mathrm{vir})$ | A 块交换用 |
| `ri_ov_exch` | $B_{ia, P}$ | $(P i, a)$ | $(n_\mathrm{aux} n_\mathrm{occ}, n_\mathrm{vir})$ | B 块交换用 |
| `fxc` | `FXCMatvecData` | — | 见下文 | XC 核表 (singlet) |

## 函数 `a_matvec`

函数路径：`ri_tddft::matvec::a_matvec`

计算 A 块矩阵-矢量积 $K^\mathbb{A}_{ia} = \sum_{jb} A_{ia,jb} z^\mathbb{A}_{jb}$，分四步：

$$
\begin{aligned}
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} (\varepsilon_a - \varepsilon_i)\, z^\mathbb{A}_{ia}
&& \text{(eq.1 对角)} \\
\mathscr{T}_{P}^{\mathbb{A}} &= \sum_{jb} B_{jb, P}\, z_{jb}^{\mathbb{A}}
&& \text{(eq.2)} \\
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} \kappa_c \sum_{P} B_{ia, P}\, \mathscr{T}_{P}^{\mathbb{A}}
&& \text{(eq.3 Coulomb)} \\
\mathscr{T}_{(Pa), j}^{\mathbb{A}} &= \sum_{b} B_{ab, P}\, z_{jb}^{\mathbb{A}}
&& \text{(eq.4)} \\
K^\mathbb{A}_{ia} &\mathrel{{+}{=}} -c_x \sum_{jP} B_{ij, P}\, \mathscr{T}_{(Pa), j}^{\mathbb{A}}
&& \text{(eq.5 交换)}
\end{aligned}
$$

其中 $\kappa_c$ 为库仑耦合因子（单重态 2、非极化 'R' 1、三重态 0）。eq.2–3 与 eq.4–5 各为一次 DGEMV 与一次 DGEMM 链；fxc 贡献由 `fxc_matvec` 补充 (见下文)。B 块交换的指标排序 $(ib|aj)$ 由同一对 `ri_ov_exch` 张量配合转置后的振幅实现，此处不赘述。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `z` | $z_{ia}^{\mathbb{A}}$ | $(i, a)$ | $(n_\mathrm{occ} n_\mathrm{vir})$ | 列优先，与振幅打印一致 |
| `ri_ov` | $B_{ia,P}$ | $(P, ia)$ | $(n_\mathrm{aux}, n_\mathrm{occ}n_\mathrm{vir})$ | |
| `ri_oo_exch` | $B_{ij,P}$ | $(jP, i)$ | $(n_\mathrm{occ}n_\mathrm{aux}, n_\mathrm{occ})$ | |
| `ri_vv_exch` | $B_{ab,P}$ | $(Pa, b)$ | $(n_\mathrm{aux}n_\mathrm{vir}, n_\mathrm{vir})$ | |

| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| fixed | (eq.2) | $\mathscr{T}_{P}^{\mathbb{A}}$ | $(P)$ | $n_\mathrm{aux}$ | |
| fixed | (eq.4) | $\mathscr{T}_{(Pa),j}^{\mathbb{A}}$ | $(Pa, j)$ | $n_\mathrm{aux} n_\mathrm{vir} n_\mathrm{occ}$ | 交换步骤的临时量 |
| fixed | (eq.5) | $K_{ia}^{\mathbb{A}}$ | $(i, a)$ | $n_\mathrm{occ} n_\mathrm{vir}$ | |

## 函数 `b_matvec`

函数路径：`ri_tddft::matvec::b_matvec`

$$
\begin{aligned}
K^\mathbb{A}_{ia} &= \kappa_c \sum_{jb} (ia|jb)\, z_{jb}^{\mathbb{A}} - c_x \sum_{jb} (ib|aj)\, z_{jb}^{\mathbb{A}} + f^{\mathrm{xc}}_{ia,jb} z_{jb}^{\mathbb{A}}
&& \text{(eq.1)}
\end{aligned}
$$

与 A 块的差异：无对角项；Coulomb 项复用 `ri_ov` 的同一缩并；交换项的指标排序不同，由 `ri_ov_exch` 配合转置振幅完成。

## fxc 核表与 `FXCMatvecData`

函数路径：`dft::num_int::prepare_fxc_data` / `dft::num_int::fxc_matvec`

MO 模式的 XC 核以「MO 轨道值 × 核表」的方式施加。`prepare_fxc_data` 在格点上求出基态密度、调用 libxc 得到二阶核，并预先将占据/虚轨道投影到格点上：

| 成员 (`FXCMatvecData`) | 意义 | 维度大小 | 其他说明 |
|--|--|--|--|
| `nvar` | 密度变量数 | | LDA 1 / GGA 4 |
| `ngrids` | 格点数 $n_\mathrm{grid}$ | | |
| `nocc` / `nvir` | 冻结芯处理后的占据/虚轨道数 | | |
| `mo_occ` | $\varphi_i(g)$ | $(n_\mathrm{occ}, n_\mathrm{grid})$ | |
| `mo_vir` | $\varphi_a(g)$ | $(n_\mathrm{vir}, n_\mathrm{grid})$ | |
| `mo_occ_grad` / `mo_vir_grad` | $\nabla\varphi(g)$ | 各 $(n_\mathrm{occ/vir}, n_\mathrm{grid})\times 3$ | 仅 GGA |
| `wfxc` | 核表 $w(g) f_{\alpha\beta}^{\mathrm{xc}}(g)$ | $n_\mathrm{grid} n_\mathrm{var}^2$ | f 连续，$g + \alpha n_\mathrm{grid} + \beta \cdot 4 n_\mathrm{grid}$ |

核施加的数学形式（以 LDA 为例）：

$$
\begin{aligned}
\rho_z^{\mathbb{A}}(g) &= \sum_{ia} z_{ia}^{\mathbb{A}} \varphi_i(g) \varphi_a(g)
&& \text{(eq.1)} \\
v^{\mathbb{A}}(g) &= w(g) f^{\mathrm{xc}}(g)\, \rho_z^{\mathbb{A}}(g)
&& \text{(eq.2)} \\
K^\mathbb{A}_{ia} &= \sum_g \varphi_i(g)\, v^{\mathbb{A}}(g)\, \varphi_a(g)
&& \text{(eq.3)}
\end{aligned}
$$

自旋通道：MO 模式的 `wfxc` 表固定为单重态核 $f_s = 2 f_u$（`SINGLET_FXC_FACTOR`）；非极化 ('R') 响应使用裸核 $f_u$（因子 1）。**三重态在 MO 模式不受支持**——`tddft_main` 会直接报错，请使用 AO 模式（见 [tddft-ao](tddft-ao.md)）。

## 求解器接口

MO 模式的 Davidson 迭代使用逐向量接口（每次施加一列试探向量）：

- `solvers::davidson::davidson_solver`（别名 `tda_davidson_solver`），闭包类型 `FnMut(&Vec<f64>) -> Vec<f64>`；
- `solvers::davidson::lr_davidson_solver`，用于全线性响应的对称化 Casida 方程。

两者是批量接口 `davidson_solver_batched` / `lr_davidson_solver_batched` 的逐列适配包装。AO 模式使用批量接口以摊销格点求值开销（见 [tddft-ao](tddft-ao.md)）。
