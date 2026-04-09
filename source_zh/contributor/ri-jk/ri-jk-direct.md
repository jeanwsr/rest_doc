# RI-JK Direct 算法


这份文档主要目的是实现内存不充足时的 RI-JK 算法。上面我们简要回顾了内存充足时需要讨论的情形，是因为一些技巧在内存不充足时仍然可以使用。

为有效利用 RI-JK 算法，我们一般需要存储完整的 Cholesky 分解 $Y_{\mu \nu, P}$ 张量。这是非常大的三维张量，即使可以利用对称性也只能降低一半存储。对于 100 原子 TZ 基组，这种存储量就会开始变得不可接受。

但既然如此，我们就必须要对电子积分进行分批、且 Cholesky 分解也会变得困难。这里将遇到的核心困难是，是否存在完美的分批方式，以及如何具体地实现分批算法。


## 分批 RI-J 算法


**函数 `generate_vj_ri_direct_with_rstsr`**


$$
\begin{align*}D_{\mathrm{tp} (\mu \nu)}^\mathbb{A} &\mathop{\tilde{\bowtie}} D_\mathrm{\mu \nu}^\mathbb{A}
\tag{eq.1} \\
g_{\mathrm{tp}(\mu \nu), \underline{P}} &\quad \text{generate}
\tag{eq.2} \\
\mathscr{T}_{\underline{P}}^{\mathbb{A}, 1} &= \sum_{\mathrm{tp}(\mu \nu)} g_{\mathrm{tp}(\mu \nu), \underline{P}} D_{\mathrm{tp}(\mu \nu)}^\mathbb{A}
\tag{eq.3} \\
\mathscr{T}_P^{\mathbb{A}, 2} &= \sum_{Q} (\mathbf{J}^{-1})_{PQ} \mathscr{T}_P^{\mathbb{A}, 1}
\tag{eq.4} \\
g_{\mathrm{tp}(\mu \nu), \underline{P}} &\quad \text{generate}
\tag{eq.5} \\
J_{\mathrm{tp}(\mu \nu)}^\mathbb{A} &= \sum_{\underline{P}} g_{\mathrm{tp}(\mu \nu), \underline{P}} \mathscr{T}_P^{\mathbb{A}, 2}
\tag{eq.6} \\
J_{\mu \nu}^\mathbb{A} &\bowtie J_{\mathrm{tp}(\mu \nu)}^\mathbb{A}
\tag{eq.7}
\end{align*}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `dms` | $D_\mathrm{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ | sym |
| `batch_size` | | | | for index $P$ |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.2), (eq.5) | $g_{\mathrm{tp}(\mu \nu), \underline{P}}$ | $(\mathrm{tp} (\mu \nu), \underline{P})$ | $(n_\mathrm{tp(basis)}, n_\mathrm{batch})$ | (eq.2) 与 (eq.5) 独立 |
| fixed | (eq.1) | $D_{\mathrm{tp} (\mu \nu)}^\mathbb{A}$ | $(\mathrm{tp} (\mu \nu), \mathbb{A})$ | $(n_\mathrm{tp(basis)}, n_\mathrm{set})$ |
| fixed | (eq.3) | $\mathscr{T}_{P}^{\mathbb{A}, 1}$ | $(P, \mathbb{A})$ | $(n_\mathrm{aux}, n_\mathrm{set})$ |
| fixed | (eq.4) | $\mathscr{T}_{P}^{\mathbb{A}, 2}$ | $(P, \mathbb{A})$ | $(n_\mathrm{aux}, n_\mathrm{set})$ | 复用 (eq.3) $\mathscr{T}_{P}^{\mathbb{A}, 1}$ |
| fixed | (eq.4) | $(\mathbf{J}^{-1})_{PQ}$ | $(P, Q)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ | LAPACK DSYSV |
| fixed | (eq.6) | $J_{\mathrm{tp}(\mu \nu)}^\mathbb{A}$ | $(\mathrm{tp} (\mu \nu), \mathbb{A})$ | $(n_\mathrm{tp(basis)}, n_\mathrm{set})$ | 复用 (eq.1) $D_{\mathrm{tp} (\mu \nu)}^\mathbb{A}$ |
| fixed | (eq.7) | $J_{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ |


整个分批 RI-J 算法，比起内存充足时的 RI-J 算法，代价是

- 每次 SCF 迭代都需要重新计算一次电子积分；这是无法避免的代价。
- (eq.2) 与 (eq.5) 需要计算两次电子积分，即电子积分的计算代价是两倍。关于这一点，我总怀疑应该有其他办法，但暂时没有想出来。 
- (eq.4) 的计算涉及到求逆；在我们实际程序中，使用了基于 LU 分解的 `rt::linalg::solve` 函数。但需要指出，这里实际上也可以用 Cholesky 分解进行计算。但考虑到这一步的计算消耗一般不是大的 (相比于电子积分的巨大 $O(N^3)$)，我认为这一步将效率优化到底的必要性不大。

## 分批 RI-K 算法：可完整储存 $g_{\mu i, P}$ 的情景 (semi-direct)


首先，RI-K 计算一般都假定使用轨道系数，通常会更高效，实现难度也相对较小。

不同于 RI-J，RI-K 的分批算法非常容易导致计算数量级的上升：或者是张量缩并部分升到低 $O(N^5)$、或者是电子积分升到低 $O(N^4)$。这通常是比较难以接受的。

因此，我们仍然需要考虑，是否可以存储 $g_{\mu i, P}$ 的中间张量。这对避免进一步提升计算复杂度而言，是至关重要的。


**函数 `generate_vk_ri_semi_direct_coeff_with_rstsr`**


$$
\begin{align*}
\tilde{C}_{\mu i}^\mathbb{A} &= C_{\mu i} \sqrt{n_i}
\tag{eq.1} \\
L_{PQ} &= \mathrm{Cholesky} (J_{PQ})
\tag{eq.2} \\
g_{\mathrm{tp}(\mu \nu), \underline{P}} &\quad \text{generate}
\tag{eq.3} \\
g_{\mu \nu, \underline{P}} &\bowtie g_{\mathrm{tp}(\mu \nu), \underline{P}}
\tag{eq.4} \\
g_{\mu i, \underline{P}}^\mathbb{A} &= \sum_\nu g_{\mu \nu, \underline{P}} \tilde{C}_{\nu i}^\mathbb{A}
\tag{eq.5} \\
Y_{\mu i, P}^\mathbb{A} &= \sum_Q (\mathbf{L}^{-1})_{QP} g_{\mu i, Q}^\mathbb{A}
\tag{eq.6} \\
K_{\mu \nu}^\mathbb{A} &= \sum_{i P} Y_{\mu i, P}^\mathbb{A} Y_{\nu i, P}^\mathbb{A}
\tag{eq.7}
\end{align*}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `mo_coeff` | $C_{\mu p}^\mathbb{A}$ | $(\mu, p, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{orb}, n_\mathrm{aux})$ |
| `mo_occ` | $n_{p}^\mathbb{A}$ | $(p, \mathbb{A})$ | $(n_\mathrm{orb}, n_\mathrm{aux})$ |
| `batch_size` | | | | for index $P$ |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.3) | $g_{\mathrm{tp}(\mu \nu), \underline{P}}$ | $\mathrm{tp}(\mu \nu), \underline{P}$ | $(n_\mathrm{tp(basis)}, n_\mathrm{batch})$ |
| fixed | (eq.1) | $\tilde{C}_{\mu i}^\mathbb{A}$ | $(\mu, i, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{occ}, n_\mathrm{set})$ |
| fixed | (eq.2) | $L_{PQ}$ | $(P, Q)$ | $(n_\mathrm{aux}, n_\mathrm{aux})$ |
| fixed | (eq.5) | $g_{\mu i, P}^\mathbb{A}$ | $(\mu, i, P)$ | $(n_\mathrm{basis}, n_\mathrm{occ}, n_\mathrm{aux})$ | 主要内存消耗 |
| fixed | (eq.7) | $K_{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ |
| thread | (eq.4) | $g_{\mu \nu, \underline{P}}$ | $(\mu, \nu)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ |


这里计算，代价是
- 每次 SCF 迭代计算时，需要重新计算一遍 3c-2e ERI $g_{\mu \nu, P}$；这是不可避免的。
- 对每个性质角标 $\mathbb{A}$ 作循环时，需要重新计算一遍 3c-2e ERI；这或许有其他办法，但我暂时想不出。意味着 RHF 并没有额外代价，但对于 UHF 或梯度性质计算时，就会有比较严重的额外电子积分计算代价。

需要留意，上述计算中，每次都对 $\mathbb{A}$ 作循环，而不是直接存储 $g_{\mu i, P}^\mathbb{A}$ 的四维张量；存储一个四维张量还是太恐怖了。

同时留意 (eq.6) 没有额外的内存消耗。这里使用的是 BLAS 的 `DTRSV` 函数 (三角矩阵求解)，因此不需要额外内存。


### 分批 RI-K 算法：不能完整储存 $g_{\mu i, P}$ 的情景 (direct)


在无法存储所有 $g_{\mu i, P}$ 的情况下，某种程度上，最好的解决方法并非使用 RI-K，而是使用 4c-2e ERI 作 integral screening。对于 4c-2e ERI，其分批算法是比较直观的。

但对于 3c-2e ERI，中间又涉及到一步 2c-2e ERI 的求逆，所有指标都至少耦合了一次、难以拆散。因此，我们必须要接受因分批算法而导致的计算复杂度上升的代价。由于电子积分是 $O(N^3)$ (尽管是代价很大的三次方算法)，我们能接受其代价上升到 $O(N^4)$。

同时，这种情形下，使用分子轨道代替密度矩阵的劣势会出现：当分批数量较大时，电子积分与轨道系数的缩并计算代价将会是接近 $O(N^5)$ 的。因此，我们这里将要采用密度矩阵的 RI-K 分批算法。

为了程序设计上的方便，这里我们将不会利用 $g_{\mu \nu, P}$ 关于 $(\mu, \nu)$ 的二重对称性，且分批将会针对 $\mu$ 指标进行。如果分批较小，计算代价会稍微多一些但不会超过两倍；分批较多时，该算法应在 RI 前提下，很难再更快了。


**函数 `generate_vk_ri_direct_dm_with_rstsr`**


$$
\begin{align*}
L_{PQ} &= \mathrm{Cholesky} (J_{PQ})
\tag{eq.1} \\
g_{\underline{\mu} \nu, P} &\quad \text{generate}
\tag{eq.2} \\
\mathscr{G}_{\underline{\mu} \nu, P}^\mathbb{A} &= \sum_\lambda g_{\underline{\mu} \nu, P} D_{\nu \lambda}^\mathbb{A}
\tag{eq.3} \\
\mathscr{T}_{\underline{\mu} \nu, P}^\mathbb{A} &= \sum_Q (\mathbf{L}^{-1})_{QP} \mathscr{G}_{\underline{\mu} \nu, Q}^\mathbb{A}
\tag{eq.4} \\
\mathscr{Y}_{\underline{\mu} \nu, P}^\mathbb{A} &= \sum_Q (\mathbf{L}^{-1})_{PQ} \mathscr{T}_{\underline{\mu} \nu, Q}^\mathbb{A}
\tag{eq.5} \\
g_{\underline{\nu} \kappa, P} &\quad \text{generate}
\tag{eq.6} \\
K_{\underline{\mu \nu}}^\mathbb{A} = K_{\underline{\nu \mu}}^\mathbb{A} &= \sum_{\kappa P} \mathscr{Y}_{\underline{\mu} \kappa, P}^\mathbb{A} g_{\underline{\nu} \kappa, P}
\tag{eq.7}
\end{align*}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `dms` | $D_\mathrm{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ | sym |
| `batch_size` | | | | for index $\mu$ of $g_{\mu \nu, P}$ |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.2) | $g_{\underline{\mu} \nu, P}$ | $(\underline{\mu}, \nu, P)$ | $(n_\mathrm{batch}, n_\mathrm{basis}, n_\mathrm{aux})$ |
| batched | (eq.5) | $\mathscr{Y}_{\underline{\mu} \nu, P}^\mathbb{A}$ | $(\underline{\mu}, \nu, P, \mathbb{A})$ | $(n_\mathrm{batch}, n_\mathrm{basis}, n_\mathrm{aux}, n_\mathrm{set})$ | (eq.3)-(eq.5) 使用同一内存 |
| batched | (eq.6) | $g_{\underline{\nu} \kappa, P}$ | $(\underline{\nu}, \kappa, P)$ | $(n_\mathrm{batch}, n_\mathrm{basis}, n_\mathrm{aux})$ | 覆盖 (eq.2) |
| fixed | (eq.1) | $L_{PQ}$ | $(P, Q)$ | $(n_\mathrm{aux}, n_\mathrm{aux})$ |
| fixed | (eq.7) | $K_{\underline{\mu \nu}}^\mathbb{A}$ | $(\underline{\mu}, \nu, \mathbb{A})$ | $(n_\mathrm{batch}, n_\mathrm{basis}, n_\mathrm{set})$ |
