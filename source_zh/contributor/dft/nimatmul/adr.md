# nimatmul 设计决定

本文记录 nimatmul 模块的主要设计决定与理由，涉及模块的 API 概念 (见 [concept.md](concept.md)) 与程序接口 (见 [basic-api.md](basic-api.md)) 背后的取舍。

:::{note}
本文是 nimatmul 模块内的设计记录，不采用 `contributor/adr` 的正式 ADR 流程 (编号、状态、决定变更管理)；正式 ADR 的撰写规范见 [ADR 使用规范与说明](../../adr/adr-rules-and-explanation.md)。本文各节先以条目块给出决定要点，随后是解释性内容。
:::

## 1. 张量库与矩阵后端 (rstsr + BLAS)

**决定**

- nimatmul 的张量运算统一使用 rstsr，矩阵运算经其 BLAS 后端 (`DeviceBLAS`) 完成。
- 不在 nimatmul 内自写矩阵乘法 kernel；热点缩并函数以标准 GEMM 表达。

**解释**

- rstsr 是 Rust 下的张量运算库 (类似 NumPy)，支持任意维度张量、列优先/行优先两种布局。使用 rstsr 而非直接操作原始数组，可以复用其完整的张量工具链。
- rstsr 的 `i()` 切片语法 (类似 NumPy 的索引) 使得多维张量的子视图操作非常简洁，这对 DFT 格点积分中频繁出现的 `[ngrids, nao, ncomp]` 三维张量的切片操作、以及更高维度的开壳层 `fxc_eff`/`kxc_eff` 尤为重要。
- rstsr 的运算与 NumPy 运算基本一一对应。这使得我们可以非常方便地将 NumPy 参考实现中的张量操作直接翻译为 rstsr 代码，降低实现难度并减少出错概率。
- 展示原型 (rstsr-showcase-dft-grids) 曾使用 Faer 后端，以避免 BLAS 依赖、方便编译。REST 统一使用 BLAS 后端 (OpenBLAS)；这带来链接与线程调度的约束：nimatmul 严格要求 openblas 链接，若同时链接 MKL 且优先级更高，算法会因线程冲突而非常慢 (`nimatmul.rs` 文件头注记)。并行 BLAS 的线程调度规范另见 [ADR-0001：并行 BLAS 线程调度](../../adr/adr-0001-blas-threads.md)。

## 2. 列优先 (column-major)

**决定**

- nimatmul 的所有张量采用列优先 (column-major) 布局。

**解释**

这是由于 REST 项目采用列优先。

实际上，大多数计算化学程序确实是采用列优先的。目前计算化学惯用语言是 Fortran/C++；其中 Fortran 默认 column-major，而 C++ 许多程序采用的 Arma 框架也是 column-major。PySCF 作为 Python 的计算化学库，虽然 Python 本身是 row-major，但 PySCF 内部的 C 代码有一些是 column-major 的；其 Python 接口，特别是 DFT 部分，经常传递的是 F-contiguous 或混合连续的 NumPy 高维张量。

## 3. 密度矩阵的对称性与 contract_ao_wv 的系数

**决定**

- `contract_ao_wv_without_symmetrize` 先计算非对称的半矩阵，各密度类型的贡献累加到同一个半矩阵中，最后一步以 $V \leftarrow V + V^T$ 统一对称化。
- 缩并系数：LDA 贡献 0.5，GGA 贡献 1.0，mGGA 的 tau 贡献 0.25。

**解释**

这一策略导致不同密度类型下系数的差异：

- LDA：系数 0.5。对称化后 $V_{\mu\nu} = 0.5 \times (\text{half}) + 0.5 \times (\text{half})^T = \text{full}$，恰好恢复完整值。
- GGA：系数 1.0。因为 GGA 贡献本身是 $\phi_{g \mu}^t (w_g f^{\rho_t} \phi_{g \nu}) + \text{swap}(\mu,\nu)$，两项恰好对应非对称半矩阵与其转置。
- mGGA (tau)：系数 0.25。对称化后得到 $0.25 \times (\text{half}) + 0.25 \times (\text{half})^T = 0.5 \times \text{full}$，与 $\frac{1}{2}\sum_t w_g f^\tau \phi_{g \mu}^t \phi_{g \nu}^t$ 一致 ($\mu, \nu$ 对称性与 $t$ 指标无关)。

这一策略的优点是：所有密度类型的贡献可以统一地累加到同一个非对称半矩阵中，最后一步对称化即可。缺点是：对于 LDA 和 mGGA，由于 $\mu, \nu$ 对称性本可减少一半计算量，当前策略实际上计算了完整的非对称矩阵再对称化，浪费了一半的 GEMM 计算量。但考虑到 GGA 是主流任务 (其 $\rho_g^t$ 项不能利用对称性)，且统一策略简化了程序实现，暂不引入 SYRK 类优化 (另见第 12 节)。

## 4. 对负占据数没有支持

**决定**

- bra-ket 形式的密度生成不支持负占据数；要求占据数 $n_i \geq 0$。

**解释**

bra-ket 形式中，bra 通常构造为 $C_{\mu i} \sqrt{n_i}$ ($n_i$ 为占据数)。这要求占据数 $n_i \geq 0$。对于大多数自洽场方法，占据数为正或零；但对于某些特殊方法 (如 fractional occupation 或某些 DFT 对稳定性分析的处理)，可能出现负占据数。负占据数出现的情形非常少，当前程序不支持负占据数。

REST 程序中，占据数以两种方式进入：生成基态密度时将 $\sqrt{n_i}$ 乘入占据轨道系数 (`make_cpks_vxc_fxc`，见 [basic-api.md](basic-api.md) 3.2 节)；闭壳层 CP-KS 响应中占据系数不乘 $\sqrt{2}$，而以输出上的 4 倍缩放系数补回 (`get_rks_response_bra`)。两者都隐含 $n_i \geq 0$ 的前提。

如果希望计算负占据数密度矩阵对应的密度格点，或者直接传入密度 (而不传入轨道)，或者作两次计算：分别得到正占据密度、负的负占据密度，然后两者相减。

## 5. 使用 $\rho_r$ 而非 $\gamma$ 作为泛函基本变量

**决定**

- 泛函计算与势缩并中，以密度梯度分量 $\rho_g^r$ ($r \in \{x,y,z\}$) 而非 $\gamma = |\nabla\rho|^2$ 作为泛函导数的基本变量。

**解释**

理由已在 [concept.md](concept.md) 2.2 节步骤二中说明：$\nabla\rho$ 是密度矩阵的一阶量，而 $\gamma$ 是二阶量。使用 $\rho_r$ 使得后续程序推导更为简单。

代价是格点维度从 LDA/GGA/mGGA 的 1/2/3 增加到 1/4/5 (自旋非极化)，但这是在 eval_xc 层面的增加，不反映在 GEMM 瓶颈中。

在程序实现中，LibXC 以 $\gamma$ 为变量；从 $\gamma$ 导数到 $\rho_r$ 导数的变换 (sigma unfolding) 在 `xceff/xc_deriv.rs` 中通过链式法则实现：

$$
\frac{\partial(f\rho)}{\partial\rho_r} = 2 f^\gamma \rho_r, \quad r \in \{x, y, z\}
$$

对于二阶和三阶导数，还需要处理 $\partial^2/\partial\gamma^2$ 的对角修正项 (`transform_xc_inner` 中对 $\chi, \chi' \in \{x, y, z\}$ 对角元加上 $2 f^{\gamma\gamma}$ 相关的贡献)。

## 6. LAPL 型 mGGA 的限制

**决定**

- 密度格点生成支持 LAPL 分量 (`XCDenType::LAPL`)，但势矩阵组装与泛函计算不支持 LAPL。

**解释**

当前程序可以计算 LAPL 密度格点：

$$
\nabla^2 \rho_g = 4 \tau_g + 2 \sum_{\mu\nu} \varphi_{g \mu} D_{\mu\nu} (\varphi_{g \nu, xx} + \varphi_{g \nu, yy} + \varphi_{g \nu, zz})
$$

但 `contract_ao_wv` 不支持 LAPL 缩并，`libxc_eval_eff` 遇到需要 Laplacian 的泛函会直接报错。这是因为 LAPL 缩并需要 AO 二阶导数 ($n_\text{comp} = 10$)，增加 GEMM 数量，且 LAPL 型泛函的应用场景有限。

## 7. XCDenType 枚举的递进式设计

**决定**

- `XCDenType` 枚举 (`RHO`/`SIGMA`/`TAU`/`LAPL`) 采用递进式设计：每个高级类型包含所有低级类型的分量。

**解释**

这使得 `XCDenType` 可以同时控制：

- 输出密度格点的分量数 `num_nvar()` (1/4/5/6)
- 所需 AO 导数阶 `num_ao_deriv()` (0/1/1/2)
- 所需 AO 分量数 `num_ao_comp()` (1/4/4/10，经 `AO_DERIV_DIM = [1, 4, 10, 20, 35]` 查表)

密度分量顺序统一为 $\rho, \rho_x, \rho_y, \rho_z, \tau, \nabla^2\rho$，与 `XCDenType` 的递进关系一致。无论泛函是否需要 $\tau$，LAPL 类型中 $\tau$ 始终是第 5 个分量 (而非第 4 个)，这保证了分量索引的一致性。

## 8. 格点分批 (nbatch) vs 格点分块 (nchunk)

**决定**

- `NIMatmul` 以 `nbatch` 控制内存分批、`nchunk` 控制并行分块，两者独立设置；`nbatch` 应为 `nchunk` 的倍数。

**解释**

`NIMatmul` 中有两个粒度参数：

- `nbatch`：内存控制参数。完整 AO 张量 `[ngrids, nao, ncomp]` 可能过大 (大体系可达数 GB)，因此按 `nbatch` 分批处理。每批独立计算 AO、密度、泛函和势矩阵，累积结果。默认为 $1536 \times 1 \times n_\text{threads}$ (线程数由 rayon 运行时决定)。
- `nchunk`：并行粒度参数。在纯函数层中，格点按 `nchunk` 分块分配给不同线程。`nchunk` 应对应 GEMM 的 KC 维度 (通常 256-512，以获得较好的缓存利用率)；默认 1536，即数倍 KC。

两者关系为 full-grid > batch > chunk > per-grid = 1。

## 9. fxc/kxc 有效势的格点空间缩并

**决定**

- 二阶及以上的泛函输出不直接进入势组装；先在格点空间与扰动密度缩并为一阶形式的“有效势”，再复用同一套 `contract_ao_wv` 缩并函数。

**解释**

对于二阶 (fxc) 和三阶 (kxc) 响应，泛函输出是高维张量 ($f^{\chi\chi'}$ 或 $f^{\chi\chi'\chi''}$)，不能直接传入 `contract_ao_wv`。需要先在格点空间作缩并：

- fxc：$\text{fxc\_eff\_contracted}^{\chi} = \sum_{\chi'} f^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}]$，得到 `[ngrids, nvar]` 的有效势
- kxc：$\text{kxc\_eff\_contracted}^{\chi} = \sum_{\chi'\chi''} f^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'] \, \xi_{g}^{\chi''}[\mathbf{R}'']$，得到 `[ngrids, nvar]` 的有效势

缩并后的有效势维度与 vxc_eff 相同，因此可以复用 `contract_ao_wv` 函数。这一设计使得 `contract_ao_wv` 的接口对所有导数阶保持一致。

## 10. UKS 的自旋维度约定

**决定**

- 开壳层张量中，自旋维度紧跟在密度变量维度之后，并置于格点维度之内侧。

**解释**

开壳层 (UKS) 的张量维度约定为：

- rho：`[ngrids, nvar, 2]` (而非 `[ngrids, 2, nvar]`)
- fxc_eff：`[ngrids, nvar, 2, nvar, 2]`
- fxc 输出：`[nao, nao, 2, nset]` (而非 `[nao, nao, nset, 2]`)，这与 PySCF 有区别。

这样约定的好处是，泛函各阶输出的维度构造是规则的 (`nvar`、`2` 逐阶重复)，且 `XCDenType` 的分量切片 (如 `rho.i((.., 1..4))`) 不受自旋维度干扰。

## 11. LDA 维度约定

**决定**

- LDA 的密度格点保留 `nvar = 1` 的维度，不作 squeeze。

**解释**

LDA 的 `nvar = 1`，即只有一个密度分量。PySCF 中经常将该分量约去 (squeeze)，使得 rho (闭壳层下) 的维度为 `[ngrids]` 而非 `[ngrids, 1]`。

但我们为了尽可能统一接口，保持所有密度类型的维度结构一致，选择不约去 LDA 的分量，使得闭壳层下 LDA 的 rho 维度为 `[ngrids, 1]`。

## 12. 性能优化的路径与暂缓项

**决定**

- 性能优化按影响范围分层进行：热点在纯函数层、数据结构在公共接口层、驱动接口保持稳定；部分优化明确暂缓。

**解释**

本模块的核心目标是 API 设计概念与正确的实现，而非性能展示；但三层架构的分离为后续优化提供了清晰的路径。以下按优化所需改动的影响范围分类讨论。

底层可直接优化的热点：

- **格点稀疏性 (non0tab)**：当前程序使用稠密 $\phi_{g \mu}$ 存储。引入 PySCF 风格的 non0tab 稀疏掩码需要在 `NIMatmul` 中增加掩码字段、在 `prepare_ao` 中生成掩码、在纯函数层增加掩码参数或新的纯函数。这影响中间层和纯函数层，但不影响驱动接口层。
- **Psi4 风格 blocking**：将 $\phi_{g \mu}$ 压缩为 $\phi_{g \mu'}^{\text{packed}}$ 加映射表，需要在中间层引入新的数据结构和纯函数层增加对应的缩并函数。同样不影响驱动接口层。
- **`contract_ao_wv` 系列函数**是计算瓶颈。当前使用标准 GEMM。考虑到 DFT 格点与基组存在稀疏性，且格点数通常远大于原子轨道数；未来可以开发专门针对格点-基组乘积的 micro-kernel，以更好地利用缓存和 SIMD 指令。对于小型体系，GEMM 已经非常高效；但对于大型体系，定制化的 kernel 可能带来显著性能提升。

需要多层协同的优化：

- **bra-ket 低秩优化**：当前已支持 `homogeneous_braket` 和 `one_bra_mult_ket` 等形式。对于 fxc 的 bra-trans 变种 (`make_fxc_pot_with_eff_bra_trans`)，输出从 `[nao, nao]` 变为 `[nao, nocc]`，这对上层算法 (如 TDDFT 求解器) 的数据结构有影响。因此，bra-trans 的推广需要驱动接口层和调用方共同适配。
- **格点分批策略**：`nbatch` 的选择影响内存用量与并行效率。对于特别大的体系，可能需要更精细的分批策略 (如按原子分批而非按格点顺序分批)，这需要在中间层调整 `NIMatmul` 的分批逻辑。
- **泛函计算并行**：`libxc_eval_eff_parallel` 的默认分块大小依密度类型调整。如果嵌套并行 (驱动接口层已在 rayon 线程池内)，需要退化为串行。未来可能需要更灵活的并行策略 (如独立线程池)。

暂不考虑的优化：

- **$\mu, \nu$ 对称性利用 (SYRK 类优化)**：$\rho_g$ 和 $\tau_g$ 的计算可以利用密度矩阵的对称性减少约一半计算量，但 $\rho_g^t$ 不行。考虑到 GGA/mGGA 是主流计算任务，且 SYRK 的 micro-kernel 与 GEMM 不同，暂不优先考虑 (见第 3 节)。
- **复数类型支持**：当前仅支持 `f64`。复数类型需要新的 micro-kernel 和共轭关系处理，暂不实现。
- **LAPL 型 mGGA 的势缩并**：见第 6 节。
