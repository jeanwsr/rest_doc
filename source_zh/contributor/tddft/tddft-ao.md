# TDDFT AO 模式：跃迁密度核

本页描述 `tddft_mode = "ao"` 下的矩阵-矢量积实现。记号沿用 [ri-jk 文档](../ri-jk/index.md)与 [tddft-mo](tddft-mo.md) 的约定：$\mu\nu$ 为原子轨道指标，$i,j$ 占据 MO，$a,b$ 虚轨道，$P$ 辅助基，$g$ 格点，$\alpha,\beta$ 密度变量分量，$\mathbb{A}$ 为试探向量集合指标，下划线角标 (如 $\underline{g}$) 表示分批指标。

## 概述

MO 模式需要预先存储三个 MO 基 RI 张量，总量为 $n_\mathrm{aux}(n_\mathrm{occ}n_\mathrm{vir} + n_\mathrm{occ}^2 + n_\mathrm{vir}^2)$；当 $n_\mathrm{occ}, n_\mathrm{vir} \gg n_\mathrm{basis}$ 时远超 AO 基的 $n_\mathrm{aux} n_\mathrm{basis}^2$（即 `scf.rimatr`），且 MO 变换本身是一次性的大 DGEMM。AO 模式的策略是：

- Davidson 迭代仍停留在 MO 振幅空间（子空间规模小）；
- 每次矩阵-矢量积，将整块试探向量变换为 **AO 跃迁密度** $P_{\mu\nu}^{\mathbb{A}}$；
- 库仑与交换项直接调用 `ri_jk` 的密度驱动接口（复用 SCF 已有的 `scf.rimatr`，无新的大张量）；
- XC 核通过 `dft::numint_matmul::NIMatmul` 在数值格点上对整块试探向量一次性批量求值（借鉴 PySCF `_gen_tda_operation` 的 `vind(zs)` 设计）；
- 结果经 MO 系数收缩回振幅空间。

## 函数 `prepare_ao_data_with_spin`

函数路径：`ri_tddft::tddft::prepare_ao_data_with_spin`（签名 `(scf, tddft_spin: Option<&str>)`；`tddft_solver` 传入输入卡解析的自旋通道，稳定性模块则无视输入卡直接传 `"singlet"`/`"triplet"`）

构造 AO 模式的 `TDDFTData`。MO 专属成员 (`fxc`, `fxc_u`, `ri_terms`) 为空/`None`，AO 侧成员如下（非限制参考下每自旋扇区一项，`Vec` 长度 = `n_sectors()`）：

| 成员 (`TDDFTData`) | 意义 | 维度大小 | 其他说明 |
|--|--|--|--|
| `c_occ` | $C_{\mu i}$ | $(n_\mathrm{basis}, n_\mathrm{occ})$ | 每扇区；空扇区 (occ = 0) 为零列矩阵 |
| `c_vir` | $C_{\mu a}$ | $(n_\mathrm{basis}, n_\mathrm{vir})$ | 每扇区 |
| `ni` | `NIMatmul` 数值积分器 | — | libcint AO 缓存；真实格点权重 |
| `fxc_eff` | 原始 XC 核（含单重态因子） | 限制性 $(n_\mathrm{grid}, n_\mathrm{var}, n_\mathrm{var})$；非限制 $(n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2)$ | 权重未乘，见下文 |
| `den_type` | `RHO` / `SIGMA` | — | 决定 $n_\mathrm{var}$ |
| `grid_batch` | 格点分批开关 | — | |
| `fxc_driver` | `Option<FxcDriver>` | — | `"mo"`/`"semitrans"`/`"dm"` 之一；`None` 仅限 HF 参考 (无核，J/K only) |
| `psi_occ` | 占据 MO 格点投影 $\psi_i(g)$ | $[n_\mathrm{grid}, n_\mathrm{occ}]$ | 仅 `MO`/`SEMITRANS` fxc 驱动构建；格点分批投影，布局 t-ready |
| `psi_occ_grad` | $\partial_d \psi_i(g)$ | $[3, n_\mathrm{grid}, n_\mathrm{occ}]$ | 仅 GGA；前导 $d$ 轴使批切片连续 |

两个早退分支：

- **HF 参考**（无 libxc 组分）：`fxc_eff`/`ni`/`fxc_driver` 全为 `None`，矩阵-矢量积只跑 RI J/K 部分。这使 AO 模式 TDDFT 与稳定性分析对 HF 参考无需数值格点即可运行；
- **RSH 泛函**：响应交换需要短程三中心积分 `scf.rimatr_sr`（由同泛函的 SCF 建立）；缺失时直接 panic 提示重跑 SCF。

**XC 核约定**：`fxc_eff` 存储的是**未乘权重的原始核**（单重态因子 $\times 2$ 已含），而 `NIMatmul` 构造时使用真实格点权重并在 `make_fxc_pot_with_eff` 内部完成权重乘法。这与 MO 模式「`wfxc` 预乘权重」的约定不同，但两者给出的最终收缩在数学上等价。

**自旋通道选择**（限制性参考由 `tddft_spin` 决定；`prepare_ao_data_with_spin` 的显式参数优先于输入卡）：

| 通道 | 核 | 求值方式 |
|--|--|--|
| 单重态 | $f_s = 2 f_u$ | 非极化求值 |
| 非极化 'R' | $f_u$ | 非极化求值（因子 1） |
| 三重态 | $f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}$ | 自旋极化求值（见下） |
| `tddft_spin = "both"` | 单重与三重各一张 | `tddft_main` 的 `run_spin` 按通道分别调用本函数，两通道核不共享 |
| 非限制参考 | 自旋分辨 $f_{\sigma_1\sigma_2}$ | 自旋极化求值于真实 $(\rho_{\alpha 0}, \rho_{\beta 0})$；`fxc_eff: [n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2]`，无单重/三重因子 |

三重态核无法由非极化求值获得：`prepare_ao_data_with_spin` 以 `LibXCSpin::Polarized` 在 $\rho_\uparrow = \rho_\downarrow = \rho/2$（GGA 时梯度取半）处求出自旋分辨核 $K[g, y_1, s_1, y_2, s_2]$（REST 的极化变换已链式法则到逐自旋梯度分量），再沿反对称方向组合：

$$
f_t[y_1, y_2](g) = \frac{1}{2} \sum_{s_1 s_2} (\pm 1)^{s_1 + s_2}\, K[g, y_1, s_1, y_2, s_2]
$$

此即 PySCF `nr_rks_fxc_st` 的三重态配方 (CPL, 256, 454)。

## 函数 `transition_density`

函数路径：`ri_tddft::matvec_ao::transition_density`

$$
P_{\mu\nu}^{\mathbb{A}} = \sum_{ia} C_{\mu i}\, z_{ia}^{\mathbb{A}}\, C_{\nu a}
$$

一次 DGEMM 完成。**注意 $P^{\mathbb{A}}$ 一般不是对称矩阵**（占据侧与虚轨道侧指标不同），这是 AO 模式所有后续算符必须面对的核心事实；XC 核与库仑项只依赖其对称部分，交换项则必须保留非对称性。

## 函数 `get_j_ao_batched`

函数路径：`ri_tddft::matvec_ao::get_j_ao_batched`（内部调用 `ri_jk::pure_incore::get_vj_ri_incore_nonsym`）

由于 $(\mu\nu|\kappa\lambda)$ 在 $(\kappa,\lambda)$ 交换下对称，库仑项只依赖跃迁密度的对称部分。折叠后走 ri-jk incore 的标准压缩缩并（记号同 [ri-jk-incore](../ri-jk/ri-jk-incore.md) eq.1–4）：

$$
\begin{aligned}
D_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} &=
\begin{cases}
P_{\mu\nu}^{\mathbb{A}} + P_{\nu\mu}^{\mathbb{A}}, & \mu \neq \nu \\
P_{\mu\mu}^{\mathbb{A}}, & \mu = \nu
\end{cases}
&& \text{(eq.1 折叠 } \tilde\bowtie\text{)} \\
\mathscr{T}_{P}^{\mathbb{A}} &= \sum_{\mathrm{tp}(\mu\nu)} D_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} Y_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.2)} \\
J_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} &= \sum_{P} \mathscr{T}_{P}^{\mathbb{A}} Y_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.3)} \\
J_{\mu\nu}^{\mathbb{A}} &\bowtie J_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}}
&& \text{(eq.4)}
\end{aligned}
$$

整块试探向量打包为 $[\mu, \nu, \mathbb{A}]$ 一次调用。与 SCF 的 `get_vj_ri_incore`（要求对称密度、折叠 $2D - \mathrm{diag}$）不同，`get_vj_ri_incore_nonsym` 的折叠以上式为准，适用于非对称跃迁密度。

## 函数 `get_k_ao_batched`

函数路径：`ri_tddft::matvec_ao::get_k_ao_batched`（内部调用 `ri_jk::pure_incore::get_vk_ri_incore_dm`）

交换项通过 `ri_jk` 的密度驱动 incore 算法计算，该算法**天然处理非对称密度**（逐辅助列 $M_{\mu\nu,\underline{P}} \cdot D \cdot M_{\nu\mu,\underline{P}}$，无折叠步骤；方程见 [ri-jk-incore](../ri-jk/ri-jk-incore.md) `get_vk_ri_incore_dm` 一节）：

$$
K_{\mu\nu}^{\mathbb{A}} = \sum_{\kappa \underline{P}} \mathscr{T}_{\mu\kappa, \underline{P}}^{\mathbb{A}}\, Y_{\nu\kappa, \underline{P}} \qquad
\mathscr{T}_{\mu\kappa, \underline{P}}^{\mathbb{A}} = \sum_{\lambda} Y_{\mu\lambda, \underline{P}} D_{\kappa\lambda}^{\mathbb{A}}
$$

**B 块**的交换指标排序 $(ib|aj)$ 由转置密度实现：$K\!\left[D^{\mathbb{A}^{\mathrm{T}}}\right]$，即交换 $c_\mathrm{occ}/c_\mathrm{vir}$ 并传入转置振幅。

**驱动方式选择** `tddft_ao_rik_driver`（仅 AO 模式）：

| 取值 | 路径 | 说明 |
|--|--|--|
| `"semitrans"` (默认) | `get_vk_ri_incore_coeff_pair` | 占据侧半转换后收缩 (见下)，**精确**且 $O(n_\mathrm{occ})$ 复杂度 |
| `"dm"` | `get_vk_ri_incore_dm` | 精确、整块批量、逐辅助列密度驱动 |
| `"lowrank"` | `get_vk_ri_incore_dm_lowrank` | 逐向量 SVD 低秩分解 (阈值 `tddft_svd_tol`)，有损近似 |

**`"semitrans"` 驱动**：振幅先折叠进占据侧系数（仅一侧变换到 MO 基，"半转换"），跃迁密度按结构精确分解 (秩 $\le n_\mathrm{occ}$，无需 SVD)：

$$
\begin{aligned}
C\!X_{\mu i}^{\mathbb{A}} &= \sum_{a} C_{\mu a}\, z_{ia}^{\mathbb{A}}
&& \text{(eq.1 一次批量 DGEMM)} \\
K_{\mu\nu}^{\mathbb{A}} &= \sum_{\underline{P},\,i} \bigl(M_{\underline{P}}\, C\!X^{\mathbb{A}}\bigr)_{\mu i}\,\bigl(M_{\underline{P}}\, C_{occ}\bigr)_{\nu i}
&& \text{(eq.2 get\_vk\_ri\_incore\_coeff\_pair)}
\end{aligned}
$$

eq.2 的右半转换 $M_{\underline{P}} C_{occ}$ 每批只算一次、跨全部试探向量复用。由于 $M_{\underline{P}}$ 对称，$K[P^{\mathrm{T}}] = K[P]^{\mathrm{T}}$：B 块直接取转置，两侧都保持 $k = n_\mathrm{occ}$。交换计算量从 $O(n_\mathrm{aux} n_\mathrm{basis}^3)$ 降为 $O(n_\mathrm{aux} n_\mathrm{basis}^2 n_\mathrm{occ})$。实测 (C6H6/PBE0 TDA, 4 线程)：DZ 70.4 → 62.9 s (−11%)，TZ 400.1 → 286.6 s (−28%) 且峰值内存省约 230 MB；能量与 `"dm"` 一致 ($\sim 10^{-14}$ Ha)。

**`"lowrank"` 驱动** `get_vk_ri_incore_dm_lowrank`（阈值 `tddft_svd_tol`，默认 1e-6）：

$$
z^{\mathbb{A}} = U \Sigma V^{\mathrm{T}} \qquad
k = \#\{\, \sigma_i \;\ge\; \varepsilon_{\mathrm{svd}}\, \sigma_{\max} \,\}
$$

交换计算量 $O(n_\mathrm{aux} n_\mathrm{basis}^2 k)$。注意这是**有损近似**：对典型体系 (如苯/def2-SVP) 试探向量矩阵接近满秩，低秩分解不带来加速；仅当跃迁密度 genuinely 低秩时有益。一般场景建议用精确的 `"semitrans"` 驱动替代。

## 函数 `fxc_matvec_ao_batched`

函数路径：`ri_tddft::matvec_ao::fxc_matvec_ao_batched`

对整块试探向量一次性完成 XC 核施加（PySCF `vind(zs)` 风格）。跃迁密度先对称化 $P^{\mathrm{sym}} = (P + P^{\mathrm{T}})/2$（`make_rho_from_dm` 的 SIGMA 响应假定对称密度；核响应在 $\mu \leftrightarrow \nu$ 对称化下不变，故为精确操作），随后：

$$
\begin{aligned}
\rho^{\mathbb{A}}_{\alpha}(g) &= \sum_{\mu\nu} P^{\mathrm{sym},\mathbb{A}}_{\mu\nu}\, \varphi_\mu(g)\, \varphi_\nu(g)
&& \text{(eq.1 make\_rho\_from\_dm)} \\
F_{\mu\nu}^{\mathbb{A}} &= \sum_{g,\alpha\beta} \varphi_\mu(g)\, w(g) f_{\alpha\beta}^{\mathrm{xc}}(g)\, \rho^{\mathbb{A}}_\alpha(g)\, \varphi_\nu(g)
&& \text{(eq.2 make\_fxc\_pot\_with\_eff)} \\
K^{\mathbb{A}}_{ia} &= \sum_{\mu\nu} C_{\mu i}\, F_{\mu\nu}^{\mathbb{A}}\, C_{\nu a}
&& \text{(eq.3 contract\_back)}
\end{aligned}
$$

eq.1 的输出为 $[n_\mathrm{grid}, n_\mathrm{var}, n_\mathrm{set}]$，eq.2 的输出为 $[n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set}]$（内部已对称化），两步都对整块 $n_\mathrm{set}$ 向量一次完成。库仑与交换是浮点量受限 (flop-bound) 的缩并，保持逐向量调用；fxc 是内存/带宽受限的格点收缩，批量化后格点 AO 值与 XC 核只读一次。

**fxc 驱动方式选择** `tddft_fxc_driver`（仅 AO 模式，默认 `"semitrans"`；`"dm"` 即上面的 eq.1–3 路径）：

| 取值 | 每次矩阵-矢量积的主导开销 | 内存开销 | 说明 |
|--|--|--|--|
| `"semitrans"` (默认) | $O(n_\mathrm{occ} n_\mathrm{basis} n_\mathrm{grid})$（振幅折入 $C_{vir}$ 后与裸 AO 收缩） | 仅 ψ_occ 表 $(1{+}3\delta_\mathrm{GGA}) n_\mathrm{occ} n_\mathrm{grid}$ | 见下 |
| `"dm"` | $O(n_\mathrm{basis}^2 n_\mathrm{grid})$（组装 $[n_\mathrm{basis},n_\mathrm{basis},m]$ 跃迁密度） | 无额外 | eq.1–3 路径 |
| `"mo"` | $O(n_\mathrm{occ} n_\mathrm{vir} n_\mathrm{grid})$ + ψ 表流量 | ψ 表 $(n_\mathrm{occ}{+}n_\mathrm{vir})(1{+}3\delta_\mathrm{GGA}) n_\mathrm{grid}$ | 见下 |

**`"semitrans"` 驱动**（默认，实现于 `fxc_mo_matvec` 的 `st` 分支）：把 $C_{vir}$ **预先折叠进振幅**，使虚轨道侧在格点上直接与裸 AO 值收缩，全程不形成 ψ_vir 表：

$$
\begin{aligned}
\tilde z^{\mathbb{A}}_{i\mu} &= \sum_{a} z^{\mathbb{A}}_{ia}\, C_{\mu a}
&& \text{(eq.S1 每次调用一次批量 DGEMM)} \\
\rho_z(g) &= \sum_{i} \psi_i(g) \sum_{\mu} \varphi_\mu(g)\, \tilde z^{\mathbb{A}}_{i\mu}
&& \text{(eq.S2 } \psi_{occ}\text{-vecdot} \times \text{裸 AO GEMM)} \\
v_{1,\alpha}(g) &= w(g)\sum_\beta f^{\mathrm{xc}}_{\alpha\beta}(g)\,\rho_\beta(g)
&& \text{(eq.S3 与其它驱动共用)} \\
E^{\mathbb{A}}_{ia} &= \sum_{\mu} C_{\mu a} \sum_g \varphi_\mu(g)\,\psi_i(g)\, v_{1,\alpha}(g)
&& \text{(eq.S4 两次 GEMM：} [\,m n_\mathrm{occ}, n_\mathrm{basis}\,] \times C_{vir}\text{)}
\end{aligned}
$$

与 `"mo"` 驱动共用同一实现框架 `fxc_mo_matvec`（仅占据侧 ψ 表缓存、格点分批 + rayon 分块），差别只在 GEMM 操作数：`"semitrans"` 的左操作数是裸 AO $[n_\mathrm{batch}, n_\mathrm{basis}]$、右操作数是折叠振幅 $[m\, n_\mathrm{occ}, n_\mathrm{basis}]$（每调用一次 eq.S1），回投时多一次 $[m\, n_\mathrm{occ}, n_\mathrm{basis}] \times C_{vir}$ GEMM；`"mo"` 的左操作数是按批投影的 ψ_vir、右操作数是原始振幅。UHF 下核收缩沿自旋分辨 XC 核 $f[g,\alpha,\sigma_1,\beta,\sigma_2]$ 展开 (eq.S3 的 $\sigma$ 双循环)。未知取值告警并回退 `"dm"`。

**`"mo"` 驱动**：即 MO 模式 fxc 算法的 AO 移植（数学上等价，仅实现不同）——`prepare_ao_data_with_spin` 以格点分批方式预先投影并缓存 occ/vir 的 MO-on-grid 表，每次矩阵-矢量积只在 occ/vir 空间收缩，与 MO 模式的 `fxc_matvec` 相同：

$$\psi_i(g) = \sum_\mu C_{\mu i}\,\varphi_\mu(g) \qquad \psi_a(g) = \sum_\mu C_{\mu a}\,\varphi_\mu(g)$$

（GGA 另缓存 $\partial_d\psi$，布局均为 $[n_\mathrm{grid}, \cdot]$）。每次矩阵-矢量积只做 occ/vir 空间的收缩（occ 侧缓存、vir 侧见下文流式说明）：

$$
\rho_0(g) = \sum_{ia} z_{ia}\,\psi_i(g)\psi_a(g), \qquad
\rho_{d+1}(g) = \sum_{ia} z_{ia}\,(\partial_d\psi_i\,\psi_a + \psi_i\,\partial_d\psi_a)(g)
$$

$$
v_{1,\alpha}(g) = w(g)\sum_\beta f^{\mathrm{xc}}_{\alpha\beta}(g)\,\rho_\beta(g), \qquad
E_{ia} = \sum_g \Lambda^\alpha_{ia}(g)\,v_{1,\alpha}(g)
$$

其中 $\Lambda^0_{ia} = \psi_i\psi_a$、$\Lambda^{d+1}_{ia} = \partial_d\psi_i\,\psi_a + \psi_i\,\partial_d\psi_a$。输出直接是 MO 振幅，无 $[n_\mathrm{basis},n_\mathrm{basis},m]$ 中间量、无 `contract_back`。

实测权衡 (C6H6/PBE0 TDA, 4 线程)：

| | dm | mo（分块流式） |
|---|---|---|
| DZ | 70.4 s / 909 MB | **46.6 s / 910 MB** |
| TZ | 400.1 s / 2383 MB | **110.7 s / 2150 MB** |

`fxc_mo_matvec` 仅缓存小的占据侧 ψ 表（`[ng, nocc]` + 梯度）；虚轨道侧投影**按格点批流式**计算（每批做一次 AO 求值 + `C_vir` 投影，批局部缓冲用后即释放），随后批内以 rayon 并行的格点分块执行（chunk ≈ 1536，`[chunk,·]` 缓冲常驻缓存；ρ 构建先收缩 z 与虚轨道侧，使 ρ GEMM 的 K = $n_\mathrm{vir}$）。能量与 `"dm"` 一致 (~2e-14 Ha)。实测两种体系下 `"mo"` 均更快且内存相当或更低；`"dm"` 保留为更简单的回退路径（无需 NIMatmul ψ 准备）。

### `grid_batch`：格点分批

完整格点的 AO 张量 $[n_\mathrm{grid}, n_\mathrm{basis}, n_\mathrm{comp}]$ 在中等体系上即可达数百 MB 乃至 GB。`grid_batch = true`（默认）时，eq.1–2 改为按格点分批执行（`NIMatmul::split_batch`，与 RKS Hessian 的 `make_hessian_setup_batched` 同一模式）：

$$
\rho^{\mathbb{A}}_{\alpha}(\underline{g}) = \sum_{\mu\nu \in \underline{g}} P^{\mathrm{sym},\mathbb{A}}_{\mu\nu}\, \varphi_\mu(\underline{g})\, \varphi_\nu(\underline{g}) \qquad
F^{\mathbb{A}} \mathrel{{+}{=}} \text{eq.2 在 } \underline{g} \text{ 上的贡献}
$$

每批的批大小为 `NIMatmul.nbatch`（默认 $1536 \times n_\mathrm{thread}$）；代价是每次矩阵-矢量积对每批重新求值 AO（libcint 调用便宜，实测约 +15% 时间换约 −40% 峰值内存，见 `tddft_ao.md` 开发笔记）。

| 内存需求类型 | 表达式 | 指标顺序 | 内存需求量 (`grid_batch=false` / `true`) | 其他说明 |
|--|--|--|--|--|
| fixed | AO 值缓存 | $(g, \mu, c)$ | $n_\mathrm{grid} n_\mathrm{basis} n_\mathrm{comp}$ / $n_\mathrm{batch} n_\mathrm{basis} n_\mathrm{comp}$ | $c$ 为密度分量数 (1/4) |
| batched | eq.1 输出 | $(\underline g, \alpha, \mathbb{A})$ | $n_\mathrm{batch} n_\mathrm{var} n_\mathrm{set}$ | |
| batched | eq.2 输出 | $(\mu, \nu, \mathbb{A})$ | $n_\mathrm{basis}^2 n_\mathrm{set}$ | |
| fixed | $f^{\mathrm{xc}}_{\alpha\beta}(g)$ XC 核 | $(g, \alpha, \beta)$ | $n_\mathrm{grid} n_\mathrm{var}^2$ | 两种模式相同 |
| fixed | `out` | $(\mu, \nu, \mathbb{A})$ | $n_\mathrm{basis}^2 n_\mathrm{set}$ | |

## 稠密小系统路径：`build_a_ao` / `build_b_ao`

函数路径：`ri_tddft::matvec_ao::build_a_ao` / `build_b_ao`

当 $n_\mathrm{occ} n_\mathrm{vir} \le 15$ 时，`tddft_main` 直接构造完整 A（或 B）矩阵做稠密对角化。AO 模式下：

$$
A_{ia,jb} = (\varepsilon_a - \varepsilon_i)\delta_{ij}\delta_{ab} + \left[\text{eq.1--3 施加于单位块 } E_{(ia),(jb)}\right]
$$

即以单位矩阵为试探向量块，经 `ao_kernel_block`（折叠 J + K (含 RSH 短程 $K_{SR}$，复用 `scf.rimatr_sr` 与相同的 K 驱动) + 批量 fxc + 收缩回 MO）一次性得到全部核贡献列，再加对角。非限制参考的稠密路径在拼接 $[\alpha;\beta]$ 振幅空间上执行。模式分派由 `tddft::build_a`/`build_b` 封装（MO 分支为逐列 `a_matvec` 循环），求解器代码不感知模式。

## 批量 Davidson 接口

AO 模式的 Davidson 迭代使用批量接口 `solvers::davidson::davidson_solver_batched` / `lr_davidson_solver_batched`（消费侧别名 `tda_davidson_solver_batched`），闭包类型：

```rust
FnMut(&MatrixFull<f64>) -> MatrixFull<f64>   // [dim, n_set] → [dim, n_set]
```

求解器把整块试探向量交给闭包，闭包内部依次执行：跃迁密度构造（一次 DGEMM 覆盖全部 $n_\mathrm{set}$ 列）→ 批量 J/K（单次 `ri_jk` 调用）→ 批量 fxc（单次 `make_rho_from_dm` + `make_fxc_pot_with_eff`）→ 收缩回 MO。并行性位于闭包内部（rayon），子空间迭代本身保持串行，避免嵌套线程池竞争。

逐向量接口 `davidson_solver` / `lr_davidson_solver`（别名 `tda_davidson_solver`、供 `ri_bse`/`scf_io` 使用）是批量核心的逐列适配包装。
