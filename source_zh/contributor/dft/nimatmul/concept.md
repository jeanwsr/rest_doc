# DFT 格点积分 API 设计概念

本文讨论 DFT 格点积分 API 设计中的一些核心概念。希望这些概念将有助于未来 DFT 格点积分的可扩展性、性能优化、API 调用友好性等方面的设计。

本文也希望以相对比较一致的文本，厘清 DFT 程序所需要的名词、概念、公式。本文的说法未必是最正统的，但这确实是设计 API 时一贯的思路。

:::{note}
本文的公式记号与张量维度遵循 [def.md](def.md) 的约定：以 $\chi$ 表示 DFT 基本参量 $\boldsymbol{\xi}$ 的分量角标；在密度表达式中以 $r$ 表示空间指标，在响应缩并中以 $t$ 表示空间指标。张量维度约定为 column-major。

本文只处理闭壳层 (RKS) 情形的公式；开壳层的维度差异参见 def.md 第 4 节与 [basic-api.md](basic-api.md)。文中涉及的函数接口清单见 [basic-api.md](basic-api.md)，设计取舍及其理由见 [adr.md](adr.md)。
:::

## 1. 通论：自洽场关键 API 设计

### 1.1 讨论前提：能量可分及其与密度矩阵的关系

计算任务的根本变量是基函数；这些基函数可能是原子轨道、平面波、实空间格点，不过原子轨道最常见。在该项目中，习惯称基函数为 AO (atomic orbital)，但这并不意味着我们只能处理原子轨道。

**能量是可以分解的**。譬如，对于 wB97X-V 泛函，其能量可以拆分为下述部分：

$$
E[\mathbf{D}] = E_\text{nuc-repl} + E_\text{kin} + E_\text{nuc} + E_\text{J} + E_\text{K} + E_\text{srK} + E_\text{xc} + E_\text{VV10}
$$

其中有一些要点是，在假定原子核不作为变量的前提下，**能量及其分项是密度矩阵 $\mathbf{D}$ 的函数**。这句话本身很简单，但隐含了下述推论或不太平常的反例：

- 密度矩阵定义为基函数系数矩阵 $C_{\mu i}$ 的函数：

    $$
    D_{\mu\nu} = \sum_i C_{\mu i} C_{\nu i}
    $$

    因此，能量及其分项也可以是基函数系数矩阵 $\mathbf{C}$ 的函数。

- 能量不能仅是系数矩阵 $\mathbf{C}$ 的函数。它必须要能写为密度矩阵 $\mathbf{D}$ 的显函数。这里的反例是各种 orbital-optimized post-SCF 方法。

- 能量与密度的关系必须严格满足 Hartree-Fock-Roothaan 方程：

    $$
    \mathbf{V}[\mathbf{D}] \mathbf{C} = \mathbf{S} \mathbf{C} \boldsymbol{\epsilon}
    $$

    其中，Fock 矩阵 $\mathbf{V}$ 是能量对密度矩阵 $\mathbf{D}$ 的梯度：

    $$
    \mathbf{V}[\mathbf{D}] = \frac{\partial E[\mathbf{D}]}{\partial \mathbf{D}}
    $$

    这里的反例是 constraint DFT (cDFT) 或 density-corrected DFT (dcDFT) 等方法。cDFT 包含约束项；dcDFT 的非自洽特性使得其 Fock 矩阵不再是能量的梯度。尽管我们这里写程序时仍然可以利用当前的 API，但在程序设计理念上，为了简化讨论，我们暂时不考虑这些方法。

- 原子结构不在考虑范围；我们只考虑电子结构 (以基函数表示的电子密度)。因此，对于原子核排斥能 $E_\text{nuc-repl}$、以及类似于 DFT-D3 的方法，我们认为它的能量项是常数。

同时，请留意，**本文的 Fock 矩阵以 $\mathbf{V} [\mathbf{D}]$ 表示，而非通常的 $\mathbf{F} [\mathbf{D}]$**。在下一小节将明确本文档记号。

### 1.2 核心问题：能量对密度矩阵导数

涉及到自洽场方法的计算问题，很多核心的技术问题，都涉及到如何计算能量对密度矩阵的导数。

我们这里需要作下述定义：

| 导数阶 | 变量 | 变量全写 | 惯用变量记号 |
|--|--|--|--|
| 1 | $\mathbf{V}$ | $\mathbf{V} [\mathbf{D}]$ | `vxc`, `v`, `fock`, `veff` |
| 2 | $\mathbf{F}$ | $\mathbf{F} [\mathbf{D}, \mathbf{R}]$ | `fxc`, `f`, `resp` |
| 3 | $\mathbf{K}$ | $\mathbf{K} [\mathbf{D}, \mathbf{R}^1, \mathbf{R}^2]$ | `kxc`, `k` |

- nimatmul 程序中，上述三个量依 DFT 泛函导数阶数命名为 `vxc`、`fxc`、`kxc`，分别对应 $\mathbf{V}$, $\mathbf{F}$, $\mathbf{K}$ 的 XC 格点积分贡献。

- 一阶导数定义

    $$
    V_{\mu \nu} = \frac{\partial E}{\partial D_{\mu \nu}}
    $$

    一阶导数在其他文献中，通常记为 $F_{\mu \nu}$。我们这里为了区分一阶、二阶导数，特意使用了不同的记号。

- 二阶导数定义

    $$
    F_{\mu \nu} = \sum_{\kappa \lambda} \frac{\partial^2 E}{\partial D_{\mu \nu} \, \partial R_{\kappa \lambda}} R_{\kappa \lambda}
    $$

    其中 $R_{\kappa \lambda}$ 是扰动密度矩阵；它在 TDDFT 中 (Casida 方程) 经常是激发态密度矩阵，在梯度问题中 (CP-HF/KS 方程) 是导数密度矩阵。

    在其他文献中，二阶导数通常记为 $\sum_{\kappa \lambda} A_{\mu \nu, \kappa \lambda} R_{\kappa \lambda}$，或 $G_{\mu \nu} [\mathbf{R}]$。

- 三阶导数定义

    $$
    K_{\mu \nu} = \sum_{\kappa \lambda} \sum_{\kappa' \lambda'} \frac{\partial^3 E}{\partial D_{\mu \nu} \, \partial R^1_{\kappa \lambda} \, \partial R^2_{\kappa' \lambda'}} R^1_{\kappa \lambda} R^2_{\kappa' \lambda'}
    $$

    其中 $R^1_{\kappa \lambda}$ 和 $R^2_{\kappa' \lambda'}$ 是两个扰动密度矩阵。

上面的所有定义，对所有可用于自洽场计算的能量分项都适用；且可以线性加和。譬如说，对于 wB97X-V 泛函，其 Fock 矩阵可以依葫芦画瓢写为：

$$
\mathbf{V} [\mathbf{D}] = \mathbf{V}_\text{nuc-repl} + \mathbf{V}_\text{kin} + \mathbf{V}_\text{nuc} + \mathbf{V}_\text{J} + \mathbf{V}_\text{K} + \mathbf{V}_\text{srK} + \mathbf{V}_\text{xc} + \mathbf{V}_\text{VV10}
$$

因此，对于任何一个能量分项而言，它都可以定义下述的程序接口 (以 Python 伪代码示意)：

```python
class EnergyComponent:
    def energy(self, dm0: np.ndarray) -> float:
    def get_fock(self, dm0: np.ndarray) -> np.ndarray:
    def get_2nd_resp(self, dm0: np.ndarray, dm1: np.ndarray) -> np.ndarray:
    def get_3rd_resp(self, dm0: np.ndarray, dm1: np.ndarray, dm2: np.ndarray) -> np.ndarray:
```

### 1.3 技术考虑：性能与可扩展性

上述的接口设计，相信是可以实现所有重要的计算化学需求。但实际使用中，我们还需要作如下的考量与适配：

- **闭壳层与开壳层的密度矩阵维度有所差异**。闭壳层维度是 `[nao, nao]`，而开壳层维度是 `[nao, nao, 2]` (col-major)。

- 密度矩阵，特别是扰动密度矩阵，**经常是多个**而非单个。因此，对于闭壳层，传入的 `dm1` 与 `dm2` 一般要允许是 3-dim 张量、或列表的 2-dim 矩阵。传入的 `dm0` 大多数情况下是单个矩阵。

- 密度矩阵通常是低秩的。它是由占据轨道系数构成的，而占据数 $n_\text{occ}$ 通常远小于基组数 $n_\text{AO}$。即使是 CP-KS 方程求解中所需要的扰动密度 (不同于自洽场密度)，它从构造上 (或者用 SVD 等方式分解上) 也可以分解为两个 $n_\text{AO} \times n_\text{occ}$ 的矩阵的乘积。因此**应该尽量利用密度矩阵的低秩结构**，以降低计算量。
    - 另一方面，对于特别大的分子，密度矩阵本身可能有足够的稀疏程度。密度矩阵零值的稀疏性与低秩结构不能 (至少难以) 同时利用；某些情况下，密度矩阵本身的稀疏性也可以利用。但我们关注的原子体系如果不超过 100 个原子，那么密度矩阵的低秩结构通常更重要。
    - 从 API 设计的角度而言，低秩性质的利用，等同于在计算函数中传入占据轨道系数矩阵：

        ```python
        def get_fock_by_occ(self, occ_coeff: np.ndarray) -> np.ndarray:
        ```

        具体的接口形式可能会有一些变化，但核心思想是，**API 设计应该允许直接传入占据轨道系数矩阵**，以利用密度矩阵的低秩结构。我们会在后面看到，在 DFT 格点积分中，我们是如何具体地利用这一点。

## 2. DFT 格点积分的计算步骤

### 2.1 与其他能量分项的比较

DFT 格点积分，与其他计算化学的能量分项贡献，有相同之处，也有不同之处。

相同点有

- 上述总论的概念，都适用于 DFT 格点积分。
- DFT 同样涉及到是否能利用占据轨道进行优化的问题。由于引入占据轨道与否，不只是会对程序实现细节产生影响，也会将 API 设计的问题变得复杂化。
- DFT 格点积分的计算复杂度相对较低，但计算量经常不算太小，特别是对于小到中等体系而言。因此，在设计 API 时，仍然需要考虑性能因素。

不同点有

- 对于 J/K 计算，一来能量是密度矩阵的二次函数 (因而 J/K 能量不存在三阶项 $\mathbf{K}$)，二来 J/K 的二阶项计算方式与一阶项非常相似。因此，尽管 J/K 计算量大、对其作近似仍然是现在电子结构程序的研究重点之一，但它的 API 设计相对简单。
- DFT 的困难在于，其能量是密度的函数 (无法 Taylor 截断到有限阶)。且越高阶，计算量越大。如果不作合理的程序设计，其程序实现难度也会很大。

### 2.2 密度生成、泛函计算、能量导数矩阵组装

DFT 格点积分的计算过程可以清晰地分为三个步骤。这三个步骤的输入输出之间有明确的依赖关系，但它们在计算特征上有显著差异。我们在这里对闭壳层 (RKS) 情形逐一说明，并区分 $\text{RHO}$ (LDA)、$\text{SIGMA}$ (GGA)、$\text{TAU}$ (mGGA) 三种密度类型 (与 def.md 相同，我们暂不考虑 $\text{LAPL}$ 型 mGGA；LAPL 的现状见 [adr.md](adr.md))。

#### 步骤一：密度格点生成 (eval_rho)

以密度矩阵 $D_{\mu\nu}$ 为输入，原子轨道格点 $\phi_{g \mu}$ 及其空间导数 $\phi_{g \mu}^r$ (其中 $r \in \{x, y, z\}$) 为中间量，生成密度格点。定义密度变量

$$
\xi_g[\mathbf{D}] := (\rho_g, \rho_g^x, \rho_g^y, \rho_g^z, \tau_g)
$$

对于闭壳层，该变量记为 $\xi_g^{\chi}$，在程序中是 `[ngrids, nvar]` 的二维数组。对于开壳层，该变量记为 $\xi_g^{\chi \sigma}$，即增加一个自旋维度，在程序中是 `[ngrids, nvar, 2]` 的三维数组。

各分量计算公式如下：

- **密度** $\rho_g$ (所有类型都需要)

    $$\rho_g = \sum_{\mu\nu} \phi_{g \mu} D_{\mu\nu} \phi_{g \nu}$$

    算法实现为一步矩阵乘法加一步数乘约化：

    $$\bar{\phi}_{g \mu} = \sum_\nu D_{\mu\nu} \phi_{g \nu} \quad \text{(GEMM)}$$

    $$\rho_g = \sum_\mu \phi_{g \mu} \bar{\phi}_{g \mu} \quad \text{(memory bounded)}$$

- **密度梯度** $\rho_g^r$ (GGA / mGGA 需要)

    $$\rho_g^r = 2 \sum_{\mu\nu} \phi_{g \mu}^r D_{\mu\nu} \phi_{g \nu}$$

    由于 $\bar{\phi}_{g \mu}$ 已经在 $\rho_g$ 的计算中得到，因此每个梯度分量只需额外一步数乘约化 $\sum_\mu \phi_{g \mu}^r \bar{\phi}_{g \mu}$，乘以系数 2。三个分量总共的额外计算量为 memory bounded。

- **动能密度** $\tau_g$ (mGGA 需要)

    $$\tau_g = \sum_{r,\mu\nu} \frac{1}{2} \phi_{g \mu}^r D_{\mu\nu} \phi_{g \nu}^r$$

    这需要对每个 $r$ 计算一步新的矩阵乘法 $\bar{\phi}_{g \mu}^{(r)} = \sum_\nu D_{\mu\nu} \phi_{g \nu}^r$ (乘以系数 1/2)，再作数乘约化 $\sum_\mu \phi_{g \mu}^r \bar{\phi}_{g \mu}^{(r)}$。三步 GEMM 的 FLOPs 约 $6 n_\text{basis}^2 n_\text{grid}$。

**低秩优化**：当密度矩阵具有低秩结构 $D_{\mu\nu} = \sum_i C_{\mu i} C_{\nu i}$ 时 ($n_\text{occ} \ll n_\text{basis}$)，可以通过占据轨道系数 $C_{\mu i}$ (bra-ket 形式) 将计算量降低：

$$\phi_{g i} = \sum_\mu \phi_{g \mu} C_{\mu i} \quad \text{(GEMM)}$$

$$\rho_g = \sum_i \phi_{g i}^{(\text{bra})} \phi_{g i}^{(\text{ket})} \quad \text{(memory bounded)}$$

对于 def2-TZVP 级别的 3-$\zeta$ 基组，$n_\text{basis} / n_\text{occ} \sim 10$，因此计算量可降低约一个数量级。

各密度类型下，密度格点生成的计算量总结如下：

| 密度类型 | 变量数 $n_\text{var}$ | AO 导数阶 | AO 分量数 $n_\text{comp}$ | GEMM 数 (DM 输入) |
|--|--|--|--|--|
| RHO (LDA) | 1 | 0 | 1 | 1 |
| SIGMA (GGA) | 4 | 1 | 4 | 1 |
| TAU (mGGA) | 5 | 1 | 4 | 4 |

密度格点生成是计算量仅次于响应矩阵生成的步骤。

#### 步骤二：泛函计算 (eval_xc)

将密度格点 $\xi_g[\mathbf{D}]$ 代入到密度泛函 $f(\rho, \nabla\rho, \tau)$ 及其各阶偏导数，得到泛函输出量。这部分的计算量是 $O(n_\text{grid})$，不是计算瓶颈。

泛函输出量的维度定义如下 (闭壳层)：

- 一阶导数 (有效势，用于 vxc)：$f_{g}^{\chi}$，维度 `[ngrids, nvar]`
- 二阶导数 (有效核，用于 fxc)：$f_{g}^{\chi \chi'}$，维度 `[ngrids, nvar, nvar]`
- 三阶导数 (有效核，用于 kxc)：$f_{g}^{\chi \chi' \chi''}$，维度 `[ngrids, nvar, nvar, nvar]`

开壳层则在每个 `nvar` 后插入自旋维度 2，例如一阶导数为 `[ngrids, nvar, 2]`。

:::{note}
**重要设计选择**：本程序使用密度梯度分量 $\rho_g^r$ 作为泛函的基本变量，而非 $\gamma = |\nabla\rho|^2$。这是因为 $\gamma$ 是密度矩阵的二阶量 (对 $\gamma$ 作密度矩阵的额外导数不为零)，而 $\nabla\rho$ 是密度矩阵严格的一阶量。这使得后续公式推导与程序实现更为简单。代价是格点维度增大 (自旋非极化 LDA/GGA/mGGA 的变量数从 1/2/3 增加到 1/4/5)；但由于 DFT 格点积分的瓶颈是 GEMM 运算而非格点维度，这一代价通常不反映在真正的计算瓶颈里。更详细的讨论见 [adr.md](adr.md)。
:::

从 $\gamma$ 到 $\rho_r$ 的变换遵循链式法则：

$$
\frac{\partial(f\rho)}{\partial\rho_r} = \frac{\partial(f\rho)}{\partial\gamma} \frac{\partial\gamma}{\partial\rho_r} = 2 f^\gamma \rho_r, \quad r \in \{x, y, z\}
$$

#### 步骤三：能量导数矩阵组装 (contract_ao_wv)

将格点数乘量 $w_g f_\text{eff}$ 与原子轨道格点 $\phi_{g \mu}$ 作缩并，得到能量导数矩阵 (Fock 矩阵的 XC 贡献)。

**一阶响应 (vxc)**。对于闭壳层情形：

$$V_{\mu\nu}^{\text{xc}} = \sum_g w_g f_g^\rho \, \phi_{g \mu} \phi_{g \nu} \quad \text{(LDA)}$$

$$V_{\mu\nu}^{\text{xc}} \leftarrow \sum_{t,g} w_g f_g^{\rho_t} \phi_{g \mu}^t \phi_{g \nu} + \text{swap}(\mu, \nu) \quad \text{(GGA)}$$

$$V_{\mu\nu}^{\text{xc}} \leftarrow \sum_{t,g} \frac{1}{2} w_g f_g^\tau \phi_{g \mu}^t \phi_{g \nu}^t \quad \text{(mGGA)}$$

上述表达式的共同结构是：**左矢 $\phi_{g \mu}^{(\text{lhs})}$、右矢 $\phi_{g \nu}^{(\text{rhs})}$、格点数乘量 $w_g f_\text{eff}$，最终对格点指标 $g$ 求和**。不同密度类型下，左矢、右矢和数乘量的具体内容不同，但结构一致。程序实现上，先计算非对称的半矩阵、再一步对称化；不同密度类型下的缩并系数有所差异，见 [adr.md](adr.md)。

**二阶响应 (fxc)**。引入 DFT 基本参量 $\boldsymbol{\xi}$ 后，fxc 的表达式为

$$
F_{\mu\nu}[\mathbf{R}] = \sum_g \sum_\chi w_g \left(\sum_{\chi'} f_g^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}]\right) \frac{\partial \xi_g^{\chi}[\mathbf{D}]}{\partial D_{\mu\nu}}
$$

其中 $\sum_{\chi'} f_g^{\chi\chi'} \xi_{g}^{\chi'}[\mathbf{R}]$ 是格点空间上的缩并，得到关于角标 $\chi$ 的格点量；随后的 $\sum_\chi (\cdots) \, \partial \xi_g^{\chi} / \partial D_{\mu\nu}$ 的处理与 vxc 完全一致。这意味着 fxc 的格点缩并步骤可以复用 vxc 的 `contract_ao_wv` 函数。

**三阶响应 (kxc)** 类似地，在格点空间缩并 $\sum_{\chi'\chi''} f_g^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'] \, \xi_{g}^{\chi''}[\mathbf{R}'']$ 后，同样复用 `contract_ao_wv`。

**计算瓶颈**。响应矩阵生成是 DFT 格点积分中计算量最大的步骤。对于 vxc，FLOPs 量级为 $O(n_\text{basis}^2 n_\text{grid} n_\text{var})$；对于 fxc，由于有 $n_\text{set}$ 个扰动密度矩阵，量级为 $O(n_\text{basis}^2 n_\text{grid} n_\text{var}^2 n_\text{set})$。密度格点生成次之。泛函计算可以忽略。

### 2.3 总体设计：分离泛函计算与格点缩并

从上一节的分析可以看出，三个计算步骤在计算特征上有显著差异：

- **密度格点生成**与**能量导数矩阵组装**都是矩阵运算 (GEMM + 数乘约化)，与泛函的具体形式无关；
- **泛函计算**是逐格点的运算 ($O(n_\text{grid})$)，与原子轨道的基组结构无关。

因此，我们将泛函计算与格点缩并 (在本程序中称为 NIMatmul) 分离为独立的模块。这一分离带来了以下好处：

1. **格点缩并函数可以接受“有效势” (eff_pot) 作为输入**，而非泛函的原始输出。有效势是泛函输出 $f_g^{\chi}$ 与密度格点 $\xi_{g}^{\chi'}[\mathbf{R}]$ (在二阶以上时) 作格点空间缩并后的格点权重向量。这使得 `contract_ao_wv` 系列函数完全不依赖 LibXC，可以独立测试与优化。

2. **泛函计算模块可以独立替换**。当前使用 LibXC，但未来可以接入 XCFun、机器学习泛函或其他泛函引擎，只要其输出格式符合有效势的约定即可。

3. **数据流清晰**。完整的数据流如下：

    ```
    // vxc
    dm → [eval_rho] → rho → [eval_xc] → vxc_eff → [contract_ao_wv] → vxc

    // fxc
    dm + dm1 → [eval_rho] → rho, rho1 → [eval_xc] → fxc_eff → [contract with rho1] → fxc_eff_contracted → [contract_ao_wv] → fxc

    // kxc
    dm + dm1 + dm2 → [eval_rho] → rho, rho1, rho2 → [eval_xc] → kxc_eff → [contract with rho1, rho2] → kxc_eff_contracted → [contract_ao_wv] → kxc
    ```

    对于一阶 (vxc)，泛函输出直接就是有效势。对于二阶 (fxc) 和三阶 (kxc)，需要在格点空间先作缩并 $\sum_{\chi'} f^{\chi\chi'} \xi_{g}^{\chi'}[\mathbf{R}]$ (或 $\sum_{\chi'\chi''} f^{\chi\chi'\chi''} \xi_{g}^{\chi'}[\mathbf{R}'] \xi_{g}^{\chi''}[\mathbf{R}'']$)，得到有效势后再进入 `contract_ao_wv`。

4. **低秩优化可以在 eval_rho 层面实现**，不影响 eval_xc 和 contract_ao_wv 的接口。本程序提供了四种密度格点生成方式 (见 [basic-api.md](basic-api.md))，以支持不同场景下的低秩优化需求。

### 2.4 三层架构：纯函数算法、DFT 公共接口、驱动接口

本程序采用三层架构设计。从底层到顶层分别是：

#### 纯函数算法层 (最底层)

这一层在 REST 中是 `numint_matmul/pure_eval_rho.rs` 与 `numint_matmul/pure_xcpot.rs` 中的函数。它们的特点是：

- **无状态**：函数不持有任何 `self` 或隐藏状态，所有输入通过参数显式传入。
- **参数简单**：输入为 tensor views 和枚举参数 (如 `XCDenType`)，输出写入预分配的 buffer (`*_with_output` 后缀)。最好不要有过于复杂的类型 (譬如完整的 grids 结构)。理想情况下，这类函数也是容易 export 到 C API 的。
- **数据结构扁平**：不涉及格点分批、AO 缓存等复杂逻辑。

纯函数也是**性能热点**所在。由于参数表完全显式，这一层可以被独立测试、替换或进一步优化，而不影响上层接口。

#### DFT 公共接口层 (中间层)

这一层是 `NIMatmul` 结构体 (`numint_matmul/nimatmul.rs`) 和泛函计算桥接 (`dft/xceff` 模块：`flags.rs`、`libxc_wrap.rs`、`xc_deriv.rs`)。它们的特点是：

- **有状态**：`NIMatmul` 管理 AO 缓存、格点分批参数 (`nchunk`、`nbatch`)、积分引擎等。
- **参数较复杂**：方法调用涉及缓存策略、分批逻辑、泛函引擎选择等。

这一层是**可能发生变化的层**。例如，更换积分引擎、改变格点存储格式 (从稠密 $\phi_{g \mu}$ 到 Psi4 blocking 或稀疏格式)、调整缓存策略等，都主要在这一层实现。但重要的是，**这一层的变动不应影响纯函数层和驱动接口层**。

#### 驱动接口层 (最上层)

这一层将三层计算步骤串联为完整的计算流程，面向最终用户。在 REST 中，这一层是 Hessian 驱动接口 (`numint_matmul/hess_rks.rs`、`hess_uks.rs`，包括 `RHessKSNIMatmul` / `UHessKSNIMatmul` 结构体及其实现的 Trait 接口)；其接口角色与 PySCF 的 `dft.numint.nr_rks` / `nr_uks` 的驱动函数类似。只要纯函数层和公共接口层的 API 不变，驱动接口层可以保持稳定。

**三层架构的核心意义**在于：性能优化应集中在纯函数层 (`contract_ao_wv` 等热点函数)，数据结构的调整、或非平凡的性能优化问题发生在中间层，而顶层用户接口保持不变。这使得不同层面的开发工作可以解耦进行。

在未来，我们也会将 DFT 层与驱动层尝试抽象出 Trait，以进一步允许各类数据结构实现相同的接口，从而增强灵活性和可扩展性。在 REST 的 Hessian 驱动层，这一设想已部分落地：`RHessKSNIMatmul` / `UHessKSNIMatmul` 通过 `HessUtilAPI`、`RHessElecInteractAPI` 等 Trait 接入 Hessian 求解器 (见 [skeleton2.md](skeleton2.md) 与 [vmat1.md](vmat1.md))。
