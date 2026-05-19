# RI-JK Incore 算法


## 内存充足时的考量因素


上述积分的表达非常简单，但需要考量的因素也有不少。


我们作分类讨论。对于内存足够的情况，我们需要考虑到

- **电子密度并非是唯一的**。上述表达式使用了单个密度矩阵的表述 ($D_{\mu \nu}$)，但实际情况是我们可能处理多密度矩阵问题：对于梯度计算问题是性质数量 $n_\mathrm{prop}$，对于自旋极化问题是待处理的自旋数量 $n_\mathrm{spin}$ (闭壳层 $n_\mathrm{spin} = 1$，非限制开壳层 $n_\mathrm{spin} = 2$，多参考与限制性开壳层可能有其他情形)。

    对于这种情况，我们需要做到的事情是，尽量避免多次电子积分的计算；如果要作分批算法，要尽量保证电子积分的生成与内存读写次数较少。

- **复数情形**。这在大多数情况下是碰不到的；但如果基组是复轨道，那么情形会有所不同。我们这里不讨论这种情况。复轨道下，上述表达式会存在指标顺序问题与共轭问题，是不正确的。

- **对称性**。电子积分、张量缩并都可以通过利用对称性得到计算量、内存占用量上的简化。这里 3c-2e ERI $g_{\mu \nu, P}$ 与密度矩阵 $D_{\mu \nu}$ 具有 $(\mu, \nu)$ 的双重对称性。

- **2c-2e ERI 的正定性**。2c-2e ERI $J_{PQ}$ 的正定性可以用于快速的电子积分分解与组合。记 Cholesky 分解后的三角矩阵 $L_{PR}$：

    $$
    J_{PQ} = \sum_R L_{PR} L_{QR}
    $$

    我们可以得到 Cholesky 分解的 3c-2e ERI $Y_{\mu \nu, P}$：

    $$
    Y_{\mu \nu, P} = \sum_{Q} g_{\mu \nu, Q} (\mathbf{L}^{-1})_{PQ}
    $$

    同时需要留意到，这里的 $(\mathbf{L}^{-1})_{PQ}$ 未必需要通过求逆计算，而是通过三角矩阵求解实现。

    但也需要留意，$J_{PQ}$ 在实际应用情景中，可能因为各种原因而表现出半正定现象。这需要用其他方式解决。

- **轨道系数替代密度矩阵**。我们可以注意到，密度矩阵可以用占据轨道表示

    $$
    D_{\mu \nu} = \sum_{i} C_{\mu i} C_{\nu i} n_i
    $$

    密度矩阵大多数时候是不满秩且半正定的 (少数例外情况包括 MP2 等微扰方法的 1-RDM)。利用这一特性，轨道系数的使用经常可以节省计算量。

    但在一些特定情况中 (特别是梯度计算时)，会遇到密度矩阵尽管不满秩，但也并非半正定的情况。例如密度矩阵的导数矩阵，是可以通过下述方式得到的：

    $$
    \mathbf{D}_{\mu \nu}^\mathbb{A} = \sum_{i} C_{\mu i} C_{\nu i}^\mathbb{A} + \text{swap} (\mu, \nu)
    $$

    其中，轨道系数导数 $C_{\nu i}^\mathbb{A}$ 是轨道系数与 U 矩阵的乘积 $\sum_{p} C_{\nu p} U_{pi}^\mathbb{A}$。对于这类型密度矩阵，在本征值的表现是正负本征值相等且各一半；在实际计算上，使用左右系数的计算方式会更自然 (参考 Psi4) 且计算量也是最省的。


## 具体算法与函数接口


- 后文中，我们记 $n_{\mathrm{tp(basis)}} = \frac{1}{2} n_\mathrm{basis} (\mathrm{basis} + 1)$，即三角存储的原子轨道矩阵大小。 
- 下划线的角标，意味着其在计算时需要被分批。注意一些算法还会对其他指标 (特别是 $\mathbb{A}$) 进行明确的循环。


**函数 `get_vj_ri_incore`**


函数路径：`ri_jk::pure_incore::get_vj_ri_incore`


$$
\begin{aligned}
D_{\mathrm{tp} (\mu \nu)}^\mathbb{A} &\mathop{\tilde{\bowtie}} D_\mathrm{\mu \nu}^\mathbb{A}
&& \text{(eq.1)} \\
\mathscr{T}_{P}^\mathbb{A} &= \sum_{\mathrm{tp} (\mu \nu)} D_{\mathrm{tp} (\mu \nu)}^\mathbb{A} Y_{\mathrm{tp} (\mu \nu), P}
&& \text{(eq.2)} \\
J_{\mathrm{tp} (\mu \nu)}^\mathbb{A} &= \sum_{P} \mathscr{T}_{P}^\mathbb{A} Y_{\mathrm{tp} (\mu \nu), P}
&& \text{(eq.3)} \\
J_{\mu \nu}^\mathbb{A} &\bowtie J_{\mathrm{tp} (\mu \nu)}^\mathbb{A}
&& \text{(eq.4)}
\end{aligned}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `cderi` | $Y_{\mathrm{tp} (\mu \nu), P}$ | $(\mathrm{tp} (\mu \nu), P)$ | $(n_\mathrm{tp(basis)}, n_\mathrm{aux})$ | $\bowtie$ |
| `dms` | $D_\mathrm{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ | sym |


**函数 `get_vk_ri_incore_coeff`**


函数路径：`ri_jk::pure_incore::get_vk_ri_incore_coeff`


$$
\begin{aligned}
\tilde{C}_{\mu i}^\mathbb{A} &= C_{\mu i} \sqrt{n_i}
&& \text{(eq.1)} \\
Y_{\mu \nu, \underline{P}} &\bowtie Y_{\mathrm{tp}(\mu \nu), \underline{P}}
&& \text{(eq.2)} \\
Y_{\mu i, \underline{P}}^\mathbb{A} &= \sum_\nu Y_{\mu \nu, \underline{P}} \tilde{C}_{\nu i}
&& \text{(eq.3)} \\
K_{\mu \nu}^\mathbb{A} &= \sum_{i \underline{P}} Y_{\mu i, \underline{P}}^\mathbb{A} Y_{\nu i, \underline{P}}^\mathbb{A}
&& \text{(eq.4)}
\end{aligned}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `cderi` | $Y_{\mathrm{tp} (\mu \nu), P}$ | $(\mathrm{tp} (\mu \nu), P)$ | $(n_\mathrm{tp(basis)}, n_\mathrm{aux})$ | $\bowtie$ |
| `mo_coeff` | $C_{\mu p}^\mathbb{A}$ | $(\mu, p, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{orb}, n_\mathrm{aux})$ |
| `mo_occ` | $n_{p}^\mathbb{A}$ | $(p, \mathbb{A})$ | $(n_\mathrm{orb}, n_\mathrm{aux})$ |
| `batch_size` | | | | for index $P$ |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.3) | $Y_{\mu i, \underline{P}}^\mathbb{A}$ | $(\mu, i, \underline{P})$ | $(n_\mathrm{basis}, n_\mathrm{occ}, n_\mathrm{batch})$ |
| fixed | (eq.1) | $\tilde{C}_{\mu i}^\mathbb{A}$ | $(\mu, i, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{occ}, n_\mathrm{set})$ |
| fixed | (eq.4) | $K_{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ |
| thread | (eq.2) | $Y_{\mu \nu, \underline{P}}$ | $(\mu, \nu)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ |


**函数 `get_vk_ri_incore_dm`**


函数路径：`ri_jk::pure_incore::get_vk_ri_incore_dm`


$$
\begin{aligned}
Y_{\mu \nu, \underline{P}} &\bowtie Y_{\mathrm{tp}(\mu \nu), \underline{P}}
&& \text{(eq.1)} \\
\mathscr{T}_{\mu \nu, \underline{P}}^\mathbb{A} &= \sum_\lambda Y_{\mu \lambda, \underline{P}} D_{\nu \lambda}^\mathbb{A}
&& \text{(eq.2)} \\
K_{\mu \nu}^\mathbb{A} &= \sum_{\kappa \underline{P}} \mathscr{T}_{\mu \kappa, \underline{P}}^\mathbb{A} Y_{\mu \kappa, \underline{P}}
&& \text{(eq.3)}
\end{aligned}
$$


| 变量名 | 变量意义 | 指标顺序 | 维度大小 | 其他说明 |
|--|--|--|--|--|
| `cderi` | $Y_{\mathrm{tp} (\mu \nu), P}$ | $(\mathrm{tp} (\mu \nu), P)$ | $(n_\mathrm{tp(basis)}, n_\mathrm{aux})$ | $\bowtie$ |
| `dms` | $D_\mathrm{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ | sym |
| `batch_size` | | | | for index $P$ |


| 内存需求类型 | 公式序号 | 表达式 | 指标顺序 | 内存需求量 | 其他说明 |
|--|--|--|--|--|--|
| batched | (eq.1) | $Y_{mu \nu, \underline{P}}$ | $(\mu, \nu, \underline{P})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{batch})$ |
| batched | (eq.2) | $\mathscr{T}_{\mu \nu, \underline{P}}^\mathbb{A}$ | $(\mu, \nu, \underline{P})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{batch})$ |
| fixed | (eq.3) | $K_{\mu \nu}^\mathbb{A}$ | $(\mu, \nu, \mathbb{A})$ | $(n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set})$ |
| thread | (eq.1) | $Y_{mu \nu, \underline{P}}$ | $(\mu, \nu)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ | 独立于其他并行区间 |
| thread | (eq.2) | $\mathscr{T}_{\mu \nu, \underline{P}}^\mathbb{A}$ | $(\mu, \nu)$ | $(n_\mathrm{basis}, n_\mathrm{basis})$ | 独立于其他并行区间 |


目前在 REST 中，已有 `SCF.generate_vj_with_ri_v_sync` 与 `SCF.generate_vk_with_ri_v` 执行相同功能，而且性能上应该不会有明显差别 (内存用量上有差别)，因此我们就不将这两个程序链入自洽场计算中了。

