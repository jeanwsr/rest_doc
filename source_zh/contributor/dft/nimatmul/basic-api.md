# nimatmul 程序接口

本文说明 nimatmul 模块的程序接口，包括三层架构 (见 [concept.md](concept.md) 2.4 节) 中纯函数算法层、DFT 公共接口层的接口清单、维度约定与关键公式，以及驱动接口层的组织方式。

本文依 REST 当前实现写成；代码位于 `rest/src/dft/numint_matmul` 与 `rest/src/dft/xceff`。

:::{note}
本文的公式记号与张量维度遵循 [def.md](def.md) 的约定。列表维度 $\mathbb{A}$ 记为 `nset`，与 def.md 第 4 节的任意性质维度 `nprop` 是同一概念：密度矩阵列表的长度，通常对应到性质计算数量，对开壳层的密度生成也经常对应到自旋 $\sigma$。

纯函数与公共接口层的密度类型参数名为 `den_type`；Hessian 驱动层 (hess_rks.rs、hess_uks.rs) 的同名参数习惯记作 `xc_type`，取值相同。

同时适用于 restricted (自旋非极化、闭壳层) 与 unrestricted (自旋极化、开壳层) 的维度，以 (R)、(U) 分别标注两种情形。
:::

## 1. 模块布局

| 文件 | 架构层 | 内容 |
|--|--|--|
| `numint_matmul/pure_eval_rho.rs` | 纯函数 | 密度格点生成 `get_rho_from_*_with_output` |
| `numint_matmul/pure_xcpot.rs` | 纯函数 | XC 势矩阵组装 `{rks,uks}_{vxc,fxc,kxc}_pot_with_eff_with_output` |
| `numint_matmul/nimatmul.rs` | 公共接口 | `NIMatmul` 结构体：AO 缓存、格点分批、密度生成与势组装方法 |
| `numint_matmul/hess_rks.rs`、`hess_uks.rs` | 驱动接口 | Hessian setup 与 CP-KS 响应；见 [skeleton2.md](skeleton2.md)、[vmat1.md](vmat1.md) |
| `dft/xceff/flags.rs` | 公共接口 | `XCDenType`、`XCSpin`、`XCPar`、`AO_DERIV_DIM` |
| `dft/xceff/libxc_wrap.rs`、`xc_deriv.rs` | 公共接口 | `libxc_eval_eff` 及 LibXC 原始输出的 $\gamma \to \rho_r$ 展开 |

:::{note}
`xceff` 模块与旧接口 `dft/libxc_itrf.rs` 的 `eval_xc_eff` 功能类似；后者服务于 `dft/num_int` 与 `dft/response` (接口栈见 [../libxc.md](../libxc.md))。两者将来可能合并重构。
:::

## 2. 纯函数算法层

### 2.1 密度格点生成 (`pure_eval_rho.rs`)

四个函数实现 [concept.md](concept.md) 2.2 节步骤一的分量公式，区别只在输入形式：密度矩阵输入作基函数指标的 GEMM，bra-ket 输入作占据轨道指标的 GEMM (低秩优化)。输出统一为 $\rho_g^{\chi \mathbb{A}}$，即 `[ngrids, nvar, nset]`。

公共参数为 `ao` (`[ngrids, nao, ncomp]` 的 AO 张量视图，`ncomp` 依 `den_type` 为 1/4/4，LAPL 为 10)、`den_type` (`XCDenType`)、`out` (预分配输出缓冲) 与 `nchunk` (并行分块大小，见 3.1 节)。

**函数 `get_rho_from_dm_with_output`**

密度矩阵输入。RHO/SIGMA 情形只需一次 GEMM (eq.1)；TAU 情形对每个空间分量 $r$ 追加一次 GEMM (eq.4 的 $\bar{\phi}_{g \mu}^{(r)}$)。**密度矩阵必须是对称的**：GEMM 只作单侧缩并，另一侧由数乘约化完成，对称性保证这与完整二次型等价。

$$
\begin{aligned}
\bar{\phi}_{g \mu} &= \sum_\nu D_{\mu \nu}^{\mathbb{A}} \, \phi_{g \nu} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_\mu \phi_{g \mu} \, \bar{\phi}_{g \mu} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= 2 \sum_\mu \phi_{g \mu}^r \, \bar{\phi}_{g \mu} && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_\mu \phi_{g \mu}^r \, \bar{\phi}_{g \mu}^{(r)} && \text{(eq.4)}
\end{aligned}
$$

其中 $\bar{\phi}_{g \mu}^{(r)} = \sum_\nu D_{\mu \nu}^{\mathbb{A}} \, \phi_{g \nu}^r$。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `dm_list` | $D_{\mu \nu}^{\mathbb{A}}$ | $(\mu, \nu)$<br>`[u, v]` | 每个 `[nao, nao]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 |
|--|--|--|--|--|
| thread | (eq.1) | $\bar{\phi}_{g \mu}$ | $(g, \mu)$ | `(nchunk, nao)` |
| thread | (eq.2)–(eq.4) | 局部输出 | $(g, \chi)$ | `(nchunk, nvar)` |

eq.1 的 $\bar{\phi}_{g \mu}$ 为 eq.2/eq.3 复用；TAU 情形逐 $r$ 重算。峰值内存为 `nthreads × nchunk × (nao + nvar)`。

**函数 `get_rho_from_homogeneous_braket_with_output`**

同构 bra-ket 输入，即 $D^{\mathbb{A}} = C^{\mathbb{A}} (C^{\mathbb{A}})^T$。占据数信息应乘入系数：闭壳层占据 2 时，用户应传入 $C_{\mu i} \sqrt{2}$ (或事后将输出乘 2；我们倾向于前者，REST 驱动层正是将 $\sqrt{n_i}$ 乘入占据轨道系数)。

$$
\begin{aligned}
\phi_{g i}^{\mathbb{A}} &= \sum_\mu \phi_{g \mu} \, C_{\mu i}^{\mathbb{A}} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_i \phi_{g i}^{\mathbb{A}} \, \phi_{g i}^{\mathbb{A}} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= 2 \sum_i \phi_{g i}^{r, \mathbb{A}} \, \phi_{g i}^{\mathbb{A}} && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_i \phi_{g i}^{r, \mathbb{A}} \, \phi_{g i}^{r, \mathbb{A}} && \text{(eq.4)}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra_list` | $C_{\mu i}^{\mathbb{A}}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- 列表中各 set 的 `nocc` 可以不同。

| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 |
|--|--|--|--|--|
| thread | (eq.1) | 轨道格点缓冲 ×2 | $(g, i)$ | `(nchunk, nocc_max)` |
| thread | (eq.2)–(eq.4) | 局部输出 | $(g, \chi)$ | `(nchunk, nvar)` |

峰值内存为 `nthreads × nchunk × (2 nocc_max + nvar)`。

**函数 `get_rho_from_one_bra_mult_ket_with_output`**

共享 bra、多 ket 的 bra-ket 输入。这是 TD/CP-KS 方程中 $X_{ai}^{\mathbb{A}}$ 类扰动的典型形式：扰动密度 $D_{\mu \nu}^{\mathbb{A}} = \sum_{ai} C_{\mu i}^{\text{bra}} X_{ai}^{\mathbb{A}} C_{\nu a}^{\text{ket}}$ 的左矢可传入占据轨道系数，右矢宜先作轨道半转换 $\tilde{C}_{\nu i}^{\mathbb{A}, \text{ket}} = \sum_a X_{ai}^{\mathbb{A}} C_{\nu a}^{\text{ket}}$ 再传入。

$$
\begin{aligned}
\phi_{g i}^{(\text{bra})} &= \sum_\mu \phi_{g \mu} \, C_{\mu i}^{(\text{bra})}, \quad \phi_{g i}^{(\text{ket}, \mathbb{A})} = \sum_\mu \phi_{g \mu} \, C_{\mu i}^{(\text{ket}, \mathbb{A})} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_i \phi_{g i}^{(\text{bra})} \phi_{g i}^{(\text{ket}, \mathbb{A})} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= \sum_i \left( \phi_{g i}^{r, (\text{bra})} \phi_{g i}^{(\text{ket}, \mathbb{A})} + \phi_{g i}^{(\text{bra})} \phi_{g i}^{r, (\text{ket}, \mathbb{A})} \right) && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_i \phi_{g i}^{r, (\text{bra})} \phi_{g i}^{r, (\text{ket}, \mathbb{A})} && \text{(eq.4)}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra` | $C_{\mu i}^{(\text{bra})}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{(\text{ket}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- `ket_list` 中每个矩阵的 `nocc` 必须与 `bra` 相同。

峰值内存为 `nthreads × nchunk × (3 nocc + nvar)`。

**函数 `get_rho_from_mult_bra_mult_ket_with_output`**

多 bra、多 ket 的 bra-ket 输入。公式与 `get_rho_from_one_bra_mult_ket_with_output` 相同，仅 bra 亦随集合 $\mathbb{A}$ 变化。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra_list` | $C_{\mu i}^{(\text{bra}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `ket_list` | $C_{\mu i}^{(\text{ket}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- `bra_list` 与 `ket_list` 长度相同，占据数一一对应。

峰值内存为 `nthreads × nchunk × (3 nocc_max + nvar)`。

### 2.2 XC 势矩阵组装 (`pure_xcpot.rs`)

所有函数共享两个私有缩并核心，即 [concept.md](concept.md) 2.2 节步骤三的 `contract_ao_wv`：

- `contract_ao_wv_without_symmetrize` 产出非对称半矩阵 (各项依 `den_type` 取舍)：

    $$
    \begin{aligned}
    (\text{half})_{\mu \nu} = \;& 0.5 \sum_g \phi_{g \mu} \, w_g f_g^{\rho} \, \phi_{g \nu} \\
    &+ \sum_{t} \sum_g \phi_{g \mu}^t \, w_g f_g^{\rho_t} \, \phi_{g \nu} \\
    &+ 0.25 \sum_{t} \sum_g \phi_{g \mu}^t \, w_g f_g^{\tau} \, \phi_{g \nu}^t
    \end{aligned}
    $$

    随后统一以 $V \leftarrow (\text{half}) + (\text{half})^T$ 对称化。系数 0.5/1.0/0.25 的理由见 [adr.md](adr.md) 第 3 节。

- `contract_ao_wv_bra` 以 bra 变换后的 AO 直接产出非对称的 `[nao, nocc]`，不作对称化，系数相应为 1.0/1.0/0.5 (见下文 bra-trans 函数的 eq.3)。

各公共函数的差别只在进入缩并核心前的“有效势”$wv$ 的构成。公共参数为 `den_type`、`ao`、`weights` ($w_g$，`[ngrids]`) 与 `nchunk`。

**函数 `rks_vxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi} &= w_g f_g^{\chi} && \text{(eq.1)} \\
V_{\mu\nu}^{\text{xc}} &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{(eq.2)}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `vxc_eff` | $f_g^{\chi}$ | $(g, \chi)$<br>`[g, x]` | `[ngrids, nvar]` |
| `vxc`<br>(output) | $V_{\mu\nu}^{\text{xc}}$ | $(\mu, \nu)$<br>`[u, v]` | `[nao, nao]` |

- `vxc_eff` 是泛函计算的输出，所用密度一般为自洽场密度。
- 输出 $V_{\mu\nu}^{\text{xc}}$ 对称。

| 内存需求类型 | 表达式 | 指标顺序 | 内存需求量 |
|--|--|--|--|
| thread | `wv` 缩并缓冲 | $(g, \mu)$ | `(nchunk, nao)` |
| thread | 局部输出 | $(\mu, \nu)$ | `(nao, nao)` |

各线程的局部输出以锁累加到全局输出。

**函数 `rks_fxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi, \mathbb{A}} &= w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_g^{\chi'}[\mathbf{R}^{\mathbb{A}}] && \text{(eq.1)} \\
F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}] &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{(eq.2)}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `fxc`<br>(output) | $F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \mathbb{A})$<br>`[u, v, A]` | `[nao, nao, nset]` |

- 输出对每个 $\mathbb{A}$ 对称。

内存需求与 `rks_vxc_pot_with_eff_with_output` 相同。

**函数 `rks_kxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi, (\mathbb{A}, \mathbb{B})} &= w_g \sum_{\chi' \chi''} f_g^{\chi\chi'\chi''} \, \xi_g^{\chi'}[\mathbf{R}'^{\mathbb{A}}] \, \xi_g^{\chi''}[\mathbf{R}''^{\mathbb{B}}] && \text{(eq.1)} \\
K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}] &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{eq.2}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \chi' \chi''}$ | $(g, \chi, \chi', \chi'')$<br>`[g, x, y, z]` | `[ngrids, nvar, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset1]` |
| `rho2` | $\xi_g^{\chi}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \mathbb{B})$<br>`[g, x, B]` | `[ngrids, nvar, nset2]` |
| `kxc`<br>(output) | $K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \mathbb{A}, \mathbb{B})$<br>`[u, v, A, B]` | `[nao, nao, nset1, nset2]` |

- `rho1` 与 `rho2` 的列表长度 `nset1`、`nset2` 可以不同。
- 输出对每个 $(\mathbb{A}, \mathbb{B})$ 对称。

**函数 `rks_fxc_pot_with_eff_bra_trans_with_output`**

bra 变换的 fxc：将 AO 先与占据轨道系数缩并，输出半转换的 `[nao, nocc]`，无需对称化。

$$
\begin{aligned}
\phi_{g i}^{(\text{bra}), *} &= \sum_\mu \phi_{g \mu}^{*} \, C_{\mu i} \quad (\text{对每个分量 } *) && \text{(eq.1)} \\
wv_g^{\chi, \mathbb{A}} &= w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_g^{\chi'}[\mathbf{R}^{\mathbb{A}}] && \text{(eq.2)} \\
F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]
&= \sum_g \phi_{g \mu} \, wv_g^{\rho, \mathbb{A}} \, \phi_{g i}^{(\text{bra})} \\
&\quad + \sum_{t,g} \left( \phi_{g \mu}^t \, wv_g^{\rho_t, \mathbb{A}} \, \phi_{g i}^{(\text{bra})} + \phi_{g \mu} \, wv_g^{\rho_t, \mathbb{A}} \, \phi_{g i}^{t, (\text{bra})} \right) \\
&\quad + \frac{1}{2} \sum_{t,g} \phi_{g \mu}^t \, wv_g^{\tau, \mathbb{A}} \, \phi_{g i}^{t, (\text{bra})} && \text{(eq.3)}
\end{aligned}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `bra` | $C_{\mu i}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `fxc`<br>(output) | $F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$<br>`[u, i, A]` | `[nao, nocc, nset]` |

- `rho1` 的 `nset` 为 ket 的数量。

| 内存需求类型 | 表达式 | 指标顺序 | 内存需求量 |
|--|--|--|--|
| persistent | $\phi_{g i}^{(\text{bra}), *}$ (eq.1) | $(g, i, *)$ | `(ngrids, nocc, ncomp)` |
| thread | `wv` 缩并缓冲 | $(g, i)$ | `(nchunk, nocc)` |
| thread | 局部输出 | $(\mu, i)$ | `(nao, nocc)` |

eq.1 的 bra 轨道格点在函数内一次性生成并全程持有。

**函数 `uks_vxc_pot_with_eff_with_output`**

公式与 `rks_vxc_pot_with_eff_with_output` 相同 (eq.1/eq.2)，$f_g^{\chi} \to f_g^{\chi \sigma}$，对每个自旋通道 $\sigma$ 独立作缩并与对称化。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `vxc_eff` | $f_g^{\chi \sigma}$ | $(g, \chi, \sigma)$<br>`[g, x, σ]` | `[ngrids, nvar, 2]` |
| `vxc`<br>(output) | $V_{\mu\nu}^{\text{xc}, \sigma}$ | $(\mu, \nu, \sigma)$<br>`[u, v, σ]` | `[nao, nao, 2]` |

**函数 `uks_fxc_pot_with_eff_with_output`**

公式与 `rks_fxc_pot_with_eff_with_output` 相同，格点空间缩并含自旋通道 $\sum_{\chi' \sigma'} f^{\chi\sigma, \chi'\sigma'} \xi^{\chi'\sigma'}[\mathbf{R}^{\mathbb{A}}]$。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \sigma \chi' \sigma'}$ | $(g, \chi, \sigma, \chi', \sigma')$<br>`[g, x, σ, y, ς]` | `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$<br>`[g, x, σ, A]` | `[ngrids, nvar, 2, nset]` |
| `fxc`<br>(output) | $F_{\mu\nu}^{\text{xc}, \sigma} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \sigma, \mathbb{A})$<br>`[u, v, σ, A]` | `[nao, nao, 2, nset]` |

- 输出对每个 $(\sigma, \mathbb{A})$ 对称。

**函数 `uks_kxc_pot_with_eff_with_output`**

公式与 `rks_kxc_pot_with_eff_with_output` 相同，各张量增加自旋维度。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \sigma \chi' \sigma' \chi'' \sigma''}$ | $(g, \chi, \sigma, \chi', \sigma', \chi'', \sigma'')$ | `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$ | `[ngrids, nvar, 2, nset1]` |
| `rho2` | $\xi_g^{\chi \sigma}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \sigma, \mathbb{B})$ | `[ngrids, nvar, 2, nset2]` |
| `kxc`<br>(output) | $K_{\mu\nu}^{\text{xc}, \sigma} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \sigma, \mathbb{A}, \mathbb{B})$ | `[nao, nao, 2, nset1, nset2]` |

**函数 `uks_fxc_pot_with_eff_bra_trans_with_output`**

公式与 `rks_fxc_pot_with_eff_bra_trans_with_output` 相同，但 bra 与输出均为 $\alpha, \beta$ 两个自旋的独立张量；由于输入与输出类型不同，不能与闭壳层函数合并 (见 3.4 节)。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \sigma \chi' \sigma'}$ | $(g, \chi, \sigma, \chi', \sigma')$ | `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$ | `[ngrids, nvar, 2, nset]` |
| `bra` | $C_{\mu i}^{(\text{bra}, \sigma)}$ | $(\mu, i)$ | 两个 `[nao, nocc_α]`、`[nao, nocc_β]` |
| `fxc`<br>(output) | $F_{\mu i}^{\text{xc}, \sigma} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$ | 两个 `[nao, nocc_σ, nset]` |

## 3. DFT 公共接口层

### 3.1 `NIMatmul` 结构体

`NIMatmul` 是格点驱动的主结构体，定义在 `nimatmul.rs` 中。其关键字段为：

| 字段 | 说明 |
|--|--|
| `cint` | 积分引擎 (libcint 包装) |
| `coords`、`weights` | 格点坐标与权重 |
| `atm_idx` | 每个格点所属的原子索引；供 Becke 格点偏置导数使用 (`usize::MAX` 表示无归属格点) |
| `quadrature_weights` | Becke 划分前的 (径向 × 角向) 积分权重；同上供格点偏置导数使用 |
| `cache_tensor` | AO 缓存，以导数阶为键 (如 `"ao_deriv0"`、`"ao_deriv1"`)，值为 copy-on-write 张量 `TsrCow` |
| `nchunk` | 并行分块大小；默认 1536，通常对应数倍 GEMM micro-kernel 的 KC 维度 |
| `nbatch` | 内存分批大小；默认 $1536 \times 1 \times n_\text{threads}$ |

`NIMatmul::new` 需要 `cint`、`coords`、`weights`、`atm_idx` 与 `quadrature_weights` 五组输入；后两组是为 Becke 格点偏置 (grid-shift) 导数准备的 (见 [becke-grid-shift.md](becke-grid-shift.md))。**格点偏置导数要求格点依原子归组**；若格点由多线程生成、原子索引交错，应先经 `regroup_grids_by_atom` 重排 (见 [becke-grid-shift.md](becke-grid-shift.md))。

**格点分批 vs 格点分块**：`nbatch` 控制内存用量 (完整 AO 张量 `[ngrids, nao, ncomp]` 可能过大)，`nchunk` 控制并行粒度。关系为 full-grid > batch > chunk > per-grid = 1。

**AO 缓存策略**：`get_cached_ao(deriv)` 在需要时计算 AO 并缓存。如果高阶导数已缓存，低阶导数可以从中切片取出，避免重复计算 (这也是缓存值采用 `TsrCow` 即 copy-on-write 类型张量的原因)。`prepare_ao(deriv)` 通过 libcint 的 `eval_gto` 计算 AO 格点，输出形状为 `[ngrids, nao, ncomp]`。`split_batch(start, end)` 从全格点实例切出一个分批实例，并同步切片已缓存的 AO 张量。

### 3.2 密度格点生成方法

| 方法 | 说明 | 常见情景 |
|--|--|--|
| `make_rho_from_dm` | 从密度矩阵列表 $D_{\mu \nu}^{\mathbb{A}}$ 生成密度 | 基础功能，post-SCF 的驰豫密度计算 |
| `make_rho_from_homogeneous_braket` | 从同构左右系数 $C_{\mu i}^{\mathbb{A}}$ 生成密度 | 自洽场 |
| `make_rho_from_one_bra_mult_ket` | 同左系数 $C_{\mu i}^{\text{bra}}$，多右系数 $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | 梯度性质 |
| `make_rho_from_mult_bra_mult_ket` | 多左系数 $C_{\mu i}^{\mathbb{A}, \text{bra}}$，多右系数 $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | |

四个方法是 2.1 节同名纯函数的封装：`&mut self` 内部完成 AO 缓存管理 (`get_cached_ao`) 与 `nchunk` 的选取，输出为新分配张量。公式、输入要求与维度见 2.1 节各函数；这里只补充 API 使用上的注意点。

**函数 `make_rho_from_dm`**

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `dm_list` | $D_{\mu \nu}^{\mathbb{A}}$ | $(\mu, \nu)$<br>`[u, v]` | 每个 `[nao, nao]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| 返回值 | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- 对于闭壳层自洽场计算，必须要将单个密度转换为长度为 1 的密度列表。
- 对于开壳层计算，可以将性质变量 $\mathbb{A}$ 视作自旋变量 $\sigma$：传入 `[dm_α, dm_β]` 两个矩阵的列表，输出 `[ngrids, nvar, 2]` 恰为开壳层密度格式。

**函数 `make_rho_from_homogeneous_braket`**

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `bra_list` | $C_{\mu i}^{\mathbb{A}}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| 返回值 | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- 不要求列表中的每个占据数是相同的。即可以传入开壳层 $\alpha, \beta$ 电子数不同的左系数矩阵。
- 闭壳层占据数通常是 2；程序要么需要用户在传入左系数 $C_{\mu i}$ 时乘以 $\sqrt{2}$，要么用户后续手动将输出乘以 2。我们倾向于希望用户采用前者的策略。

**函数 `make_rho_from_one_bra_mult_ket`**

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `bra` | $C_{\mu i}^{\text{bra}}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| 返回值 | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- 这一函数的设计初衷是为了性质计算中，在求解 TD/CP-KS 方程时，涉及到的 $X_{ai}^{\mathbb{A}}$ 在代入到 DFT 计算中时需要先作原子轨道转换；占据轨道数通常远小于原子轨道数，因此传入 $C_{\mu i}^{\text{bra}}$ 可节省计算量，另一边作轨道半转换 (见 2.1 节同名纯函数)。
- 对于开壳层情况，该函数不能同时处理 $\alpha$ 与 $\beta$ 自旋。需要分两次计算。

**函数 `make_rho_from_mult_bra_mult_ket`**

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `bra_list` | $C_{\mu i}^{\mathbb{A}, \text{bra}}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `ket_list` | $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | $(\mu, i)$<br>`[u, i]` | 每个 `[nao, nocc]`，共 `nset` 个 |
| `den_type` | | | `XCDenType` |
| 返回值 | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- `bra_list` 与 `ket_list` 长度需要相同，占据数需要一一对应。

### 3.3 泛函计算桥接 (`xceff`)

`dft/xceff/flags.rs` 定义三个枚举与一个常数：

| 名称 | 取值 | 说明 |
|--|--|--|
| `XCDenType` | `RHO` / `SIGMA` / `TAU` / `LAPL` | 密度类型；递进式设计，见 [adr.md](adr.md) 第 7 节 |
| `XCSpin` | `Unpolarized` / `Polarized` | 自旋极化与否 |
| `XCPar` | `Par { chunk_size }` / `Serial` | 泛函计算并行策略；可从 `usize`、`Option<usize>`、`bool` 转换 |
| `AO_DERIV_DIM` | `[1, 4, 10, 20, 35]` | 各 AO 导数阶的分量数 (至四阶) |

`determine_den_type` 依泛函家族决定密度类型：LDA/HybLDA → `RHO`，GGA/HybGGA → `SIGMA`，mGGA/HybMGGA → 依 `needs_laplacian()` 决定 `TAU` 或 `LAPL`；`determine_den_type_from_list` 取泛函列表中最严格的类型。

**函数 `libxc_eval_eff`**

将 LibXC 的原始输出转换为有效势格式；LibXC 以 $\gamma$ 为 GGA 变量，$\gamma \to \rho_r$ 的链式法则展开 (含二阶以上的对角修正项) 在函数内部经 `xc_deriv.rs` 的 `transform_xc_inner` 完成 (见 [adr.md](adr.md) 第 5 节)。

```rust
pub fn libxc_eval_eff(
    xc_func: &LibXCFunctional,
    rho: TsrView,
    deriv: usize,
    par: impl Into<XCPar>,
) -> Vec<Tsr>
```

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `xc_func` | | | `LibXCFunctional` |
| `rho` | (R) $\xi_g^{\chi}$ <br> (U) $\xi_g^{\chi \sigma}$ | $(g, \chi)$<br>`[g, x]` | (R) `[ngrids, nvar]` <br> (U) `[ngrids, nvar, 2]` |
| `deriv` | | | |
| `par` | | | `XCPar` |
| 返回值 | $f_g$ 及各阶导数 | | `Vec<Tsr>` |

- 生成 `LibXCFunctional` 实例时须同时指定其自旋；程序以该实例判断开闭壳层。
- `deriv` 是需要计算的**最高**导数阶；返回值含 0 至 `deriv` 所有阶的输出，第 $k$ 个元素为 $k$ 阶输出。

输出的有效势维度为：

| `deriv` | 闭壳层 | 开壳层 |
|--|--|--|
| 0 | `[ngrids]` (exc，$f_g$) | `[ngrids]` |
| 1 | `[ngrids, nvar]` (vxc_eff，$f_g^{\chi}$) | `[ngrids, nvar, 2]` |
| 2 | `[ngrids, nvar, nvar]` (fxc_eff，$f_g^{\chi \chi'}$) | `[ngrids, nvar, 2, nvar, 2]` |
| 3 | `[ngrids, nvar, nvar, nvar]` (kxc_eff，$f_g^{\chi \chi' \chi''}$) | `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |

并行分块大小的默认值依密度类型与自旋而定：RHO/非极化 16384，RHO/极化 6144，SIGMA 384，TAU/LAPL 256。**嵌套并行会退化为串行**：若调用线程已在 rayon 线程池内，泛函计算只使用单线程，以避免线程过订阅。

### 3.4 XC 势矩阵组装方法

| 方法 | 说明 |
|--|--|
| `make_vxc_pot_with_eff` | 一阶 XC 势 |
| `make_fxc_pot_with_eff` | 二阶 XC 核 |
| `make_kxc_pot_with_eff` | 三阶 XC 核 |
| `make_rks_fxc_pot_with_eff_bra_trans` | 二阶 XC 核 (bra-transformed，低秩优化) |
| `make_uks_fxc_pot_with_eff_bra_trans` | 二阶 XC 核 (bra-transformed，低秩优化，开壳层) |

五个方法是 2.2 节同名纯函数的封装：`&mut self` 内部完成 AO 缓存管理与格点权重 $w_g$ 的乘入，因此 $w_g$ 与 AO 不作为参数传入。公式与输入要求见 2.2 节各函数；这里给出统一公式与 API 使用上的注意点。

**函数 `make_vxc_pot_with_eff`**

$$
V_{\mu \nu}^\text{xc} = \sum_g \sum_\chi w_g f_g^\chi \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `vxc_eff` | (R) $f_g^{\chi}$ <br> (U) $f_g^{\chi \sigma}$ | $(g, \chi)$<br>`[g, x]` | (R) `[ngrids, nvar]` <br> (U) `[ngrids, nvar, 2]` |
| `den_type` | | | `XCDenType` |
| `spin` | | | `XCSpin` |
| 返回值 | $V_{\mu\nu}^{\text{xc}}$ | $(\mu, \nu)$<br>`[u, v]` | (R) `[nao, nao]` <br> (U) `[nao, nao, 2]` |

- `spin` 决定调用 RKS 或 UKS 纯函数。

$\partial \xi_g^{\chi} / \partial D_{\mu\nu}$ 是密度矩阵无关量，只与轨道格点 $\phi_{g \mu}$ 及其梯度有关；这与 $w_g$ 一样由 `&mut self` 管理。

**函数 `make_fxc_pot_with_eff`**

$$
F_{\mu\nu}^\text{xc} [\mathbf{R}^{\mathbb{A}}] = \sum_g \sum_\chi \left( w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}^{\mathbb{A}}] \right) \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | (R) `[ngrids, nvar, nvar]` <br> (U) `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | (R) `[ngrids, nvar, nset]` <br> (U) `[ngrids, nvar, 2, nset]` |
| `den_type`、`spin` | | | |
| 返回值 | $F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \mathbb{A})$<br>`[u, v, A]` | (R) `[nao, nao, nset]` <br> (U) `[nao, nao, 2, nset]` |

**函数 `make_kxc_pot_with_eff`**

$$
K_{\mu\nu}^\text{xc} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}] = \sum_g \sum_\chi \left( w_g \sum_{\chi' \chi''} f_g^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'^{\mathbb{A}}] \, \xi_{g}^{\chi''}[\mathbf{R}''^{\mathbb{B}}] \right) \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \chi' \chi''}$ | $(g, \chi, \chi', \chi'')$<br>`[g, x, y, z]` | (R) `[ngrids, nvar, nvar, nvar]` <br> (U) `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | (R) `[ngrids, nvar, nset1]` <br> (U) `[ngrids, nvar, 2, nset1]` |
| `rho2` | $\xi_g^{\chi}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \mathbb{B})$<br>`[g, x, B]` | (R) `[ngrids, nvar, nset2]` <br> (U) `[ngrids, nvar, 2, nset2]` |
| `den_type`、`spin` | | | |
| 返回值 | $K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \mathbb{A}, \mathbb{B})$<br>`[u, v, A, B]` | (R) `[nao, nao, nset1, nset2]` <br> (U) `[nao, nao, 2, nset1, nset2]` |

**函数 `make_rks_fxc_pot_with_eff_bra_trans`**

该函数实际上就是计算了

$$
F_{\mu i}^\text{xc} [\mathbf{R}^{\mathbb{A}}] = \sum_\nu F_{\mu\nu}^\text{xc} [\mathbf{R}^{\mathbb{A}}] C_{\nu i}
$$

但一般来说，性能最优的实现模式是对 $\partial \xi_g^{\chi} / \partial D_{\mu\nu}$ 与 $C_{\nu i}$ 的乘积进行低秩优化 (2.2 节的 `contract_ao_wv_bra`)。

之所以需要设计这样的函数，是因为 TD/CP-KS 方程中经常出现下述形式的计算问题：

$$
A_{a i, b j}^\text{xc} R_{b j} = \sum_{\mu \nu} C_{\mu a} C_{\nu i} F_{\mu\nu}^\text{xc} [\mathbf{R}]
$$

而先缩并占据轨道得到 $F_{\mu i}^\text{xc} [\mathbf{R}]$ 是性能更好的做法，因此上式化为

$$
A_{a i, b j}^\text{xc} R_{b j} = \sum_{\mu i} C_{\mu a} F_{\mu i}^\text{xc} [\mathbf{R}]
$$

得到 $F_{\mu i}^\text{xc} [\mathbf{R}]$ 这一步只是半转换。需要留意，全转换在 $n_\mathrm{set} = 1$ 的计算代价基本不会减少，在 $n_\mathrm{set} > 1$ 的情形时计算量会上升；因此对于要处理 CP-KS 的情况 (有多个性质矩阵要计算)，全转换不如半转换。因此我们在公共接口层仅提供了半转换的函数。

| 变量名 | 变量意义 | 指标顺序 | 维度大小 |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `bra` | $C_{\mu i}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| 返回值 | $F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$<br>`[u, i, A]` | `[nao, nocc, nset]` |

- `rho1` 的 `nset` 是 bra-ket 形式中 ket 的数量。
- **该函数仅用于闭壳层**。这是因为开壳层的输入 `bra` 与输出不能是单个更高一阶的、带有自旋维度的张量，而必须是两个单独的 $\alpha$ 与 $\beta$ 张量。由于类型不同，因此在 Rust 等强类型语言中无法通过一个函数实现。

**函数 `make_uks_fxc_pot_with_eff_bra_trans`**

函数意义同上，但仅适用于开壳层。`bra` 为两个自旋的左系数 `[nao, nocc_α]` 与 `[nao, nocc_β]`，返回值为两个自旋的 `[nao, nocc_σ, nset]`。由于输入与输出类型不同，因而不能与闭壳层的函数合并。

## 4. 驱动接口层

REST 中，SCF 与梯度的 XC 计算目前由 `dft/num_int` (配合 `libxc_itrf`) 承担；nimatmul 的驱动接口层是 **DFT Hessian 及其 CP-KS 响应**，位于 `hess_rks.rs` 与 `hess_uks.rs`。公式推导与实现策略见 [skeleton2.md](skeleton2.md) (能量二阶 Skeleton 导数) 与 [vmat1.md](vmat1.md) (Fock 矩阵一阶 Skeleton 导数)；这里只列举接口组织。

| 函数 / 结构体 | 说明 |
|--|--|
| `get_rho_exc_vxc_fxc` (`_uks`) | 从 AO 与 (含占据的) AO-缩并密度矩阵，得到 `rho`、`exc`、`vxc`、`fxc` |
| `eval_vxc_fxc_from_rho` (`_uks`) | 从已有的 `rho` 直接计算 `vxc_eff`、`fxc_eff` (跳过密度生成) |
| `make_cpks_vxc_fxc` (`_uks`) | 精简版 vxc/fxc 评估：bra-ket 生成密度、只算最低需要的 AO 导数阶，供 CP-KS 专用格点使用 |
| `make_hessian_setup_becke` (`_uks`) | Hessian setup 主体：格点分批循环内完成密度生成、泛函计算与各 Skeleton 中间量 |
| `get_rks_response_bra` / `get_rks_response_bra_batched` (`_uks`) | CP-KS 响应：`make_rho_from_one_bra_mult_ket` + `make_fxc_pot_with_eff_bra_trans` |
| `RHessKSNIMatmul` / `UHessKSNIMatmul` | 驱动结构体；持有 `NIMatmul`、泛函列表与中间量，通过 `HessUtilAPI`、`RHessElecInteractAPI` 等 Trait 接入 Hessian 求解器 |

驱动层的通用算法模式是对格点按 `nbatch` 分批，每批内部串联三个步骤 (以 CP-KS 响应为例，`get_rks_response_bra_batched`)：

```rust
for start in (0..ngrids).step_by(nbatch) {
    let end = (start + nbatch).min(ngrids);
    // 1. 切出分批 NIMatmul (同步切片已缓存的 AO 张量)
    let mut ni_batch = ni.split_batch(start, end);
    // 2. 密度格点生成
    let rho1 = ni_batch.make_rho_from_one_bra_mult_ket(mocc, &mo1_bra_list, den_type);
    // 3. XC 势矩阵组装 (fxc_eff 由泛函计算得到，按批切片传入)
    let resp_batch = ni_batch.make_rks_fxc_pot_with_eff_bra_trans(
        fxc_eff.i(start..end), rho1, mocc, den_type);
    resp += resp_batch;
}
```

泛函计算 (`libxc_eval_eff`) 在分批循环之外的全格点上完成一次；分批循环内只作 AO 切片、密度生成与势组装。对于 vxc 一类的非响应计算，能量与电子数的累积 ($E^\text{xc} = \sum_g w_g f_g \rho_g$、$N_e = \sum_g w_g \rho_g$) 也发生在该层。

:::{note}
闭壳层 CP-KS 响应中，占据轨道系数未乘 $\sqrt{2}$，占据数的贡献以输出上的 4 倍缩放系数补回 (`get_rks_response_bra` 中 `4.0 * resp`)。这与 3.2 节 make_rho_from_homogeneous_braket 的 $\sqrt{n_i}$ 策略是同一件事的两种实现方式。
:::
