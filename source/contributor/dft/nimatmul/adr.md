# nimatmul Design Decisions

This document records the main design decisions of the nimatmul module and their rationale, covering the API concepts of the module (see [concept.md](concept.md)) and the trade-offs behind the program interface (see [basic-api.md](basic-api.md)).

:::{note}
This document is a design record within the nimatmul section; it does not follow the formal ADR process of `contributor/adr` (numbering, status, decision-change management). For the writing conventions of formal ADRs, see `adr-rules-and-explanation.md` under `contributor/adr` (currently Chinese only). Each section below first states the decision in bullet form, followed by explanatory content.
:::

## 1. Tensor Library and Matrix Backend (rstsr + BLAS)

**Decision**

- All tensor operations of nimatmul use rstsr uniformly; matrix operations go through its BLAS backend (`DeviceBLAS`).
- No matrix-multiplication kernels are hand-written inside nimatmul; the hotspot contraction functions are expressed in terms of standard GEMM.

**Rationale**

- rstsr is a tensor-operation library for Rust (similar to NumPy), supporting arbitrary-dimensional tensors with both column-major and row-major layouts. Using rstsr rather than manipulating raw arrays directly reuses its full tensor toolchain.
- The `i()` slicing syntax of rstsr (similar to NumPy indexing) makes subview operations on multi-dimensional tensors very concise. This is especially important for the slicing of the `[ngrids, nao, ncomp]` 3-dim tensors that occur frequently in DFT grid integration, as well as for the higher-dimensional open-shell `fxc_eff`/`kxc_eff`.
- The operations of rstsr correspond essentially one-to-one with NumPy operations. This allows tensor operations from NumPy reference implementations to be translated directly into rstsr code, lowering the implementation difficulty and the chance of errors.
- The showcase prototype (rstsr-showcase-dft-grids) once used the Faer backend to avoid BLAS dependencies and ease compilation. REST uses the BLAS backend (OpenBLAS) uniformly; this brings linkage and thread-scheduling constraints: nimatmul strictly requires openblas linkage, and if MKL is linked as well with higher precedence, the algorithm becomes very slow due to thread conflicts (header note of `nimatmul.rs`). The thread-scheduling rules for parallel BLAS are given in `adr-0001-blas-threads.md` under `contributor/adr` (currently Chinese only).

## 2. Column-major

**Decision**

- All tensors of nimatmul use the column-major layout.

**Rationale**

This is because the REST project adopts column-major.

In fact, most computational-chemistry programs are indeed column-major. The customary languages of computational chemistry today are Fortran/C++; Fortran defaults to column-major, and Arma, the framework adopted by many C++ programs, is also column-major. PySCF, as a Python computational-chemistry library, runs on row-major Python; but some of the C code inside PySCF is column-major, and its Python interface — the DFT part in particular — often passes F-contiguous or mixedly-contiguous NumPy high-dimensional tensors.

## 3. Density-Matrix Symmetry and the contract_ao_wv Coefficients

**Decision**

- `contract_ao_wv_without_symmetrize` first computes the asymmetric half matrix, with the contributions of all density types accumulated into the same half matrix; a final step performs the uniform symmetrization $V \leftarrow V + V^T$.
- Contraction coefficients: LDA contribution 0.5, GGA contribution 1.0, mGGA tau contribution 0.25.

**Rationale**

This strategy leads to the coefficient differences across density types:

- LDA: coefficient 0.5. After symmetrization, $V_{\mu\nu} = 0.5 \times (\text{half}) + 0.5 \times (\text{half})^T = \text{full}$, exactly restoring the complete value.
- GGA: coefficient 1.0. The GGA contribution itself is $\phi_{g \mu}^t (w_g f^{\rho_t} \phi_{g \nu}) + \text{swap}(\mu,\nu)$; its two terms correspond exactly to the asymmetric half matrix and its transpose.
- mGGA (tau): coefficient 0.25. After symmetrization, $0.25 \times (\text{half}) + 0.25 \times (\text{half})^T = 0.5 \times \text{full}$, consistent with $\frac{1}{2}\sum_t w_g f^\tau \phi_{g \mu}^t \phi_{g \nu}^t$ (the $\mu, \nu$ symmetry is independent of the $t$ index).

The advantage of this strategy is that the contributions of all density types can be accumulated uniformly into the same asymmetric half matrix, with one final symmetrization. The disadvantage is that, for LDA and mGGA, the $\mu, \nu$ symmetry could have halved the cost, but the current strategy computes the full asymmetric matrix before symmetrizing, wasting half of the GEMM work. Considering that GGA is the mainstream workload (its $\rho_g^t$ terms cannot exploit the symmetry) and that the uniform strategy simplifies the implementation, SYRK-type optimizations are not introduced for now (see also Section 12).

## 4. No Support for Negative Occupations

**Decision**

- Bra-ket density generation does not support negative occupations; occupations $n_i \geq 0$ are required.

**Rationale**

In the bra-ket form, the bra is usually constructed as $C_{\mu i} \sqrt{n_i}$ ($n_i$ being the occupation). This requires $n_i \geq 0$. For most self-consistent-field methods the occupations are positive or zero; but for some special methods (such as fractional occupation, or certain DFT stability analyses), negative occupations can arise. Such cases are rare, and the current program does not support negative occupations.

In the REST program, occupations enter in two ways: when generating the ground-state density, $\sqrt{n_i}$ is folded into the occupied-orbital coefficients (`make_cpks_vxc_fxc`, see Section 3.2 of [basic-api.md](basic-api.md)); in the closed-shell CP-KS response the occupation coefficients are not multiplied by $\sqrt{2}$, and the contribution is restored by a factor of 4 on the output (`get_rks_response_bra`). Both implicitly assume $n_i \geq 0$.

If the density on grids corresponding to a negative-occupation density matrix is desired, either pass the density directly (instead of orbitals), or perform the computation twice: obtain the positive-occupation density and the negative-of-the-negative-occupation density separately, then subtract.

## 5. Using $\rho_r$ rather than $\gamma$ as the Fundamental Functional Variables

**Decision**

- In functional evaluation and potential contraction, the density-gradient components $\rho_g^r$ ($r \in \{x,y,z\}$) are used as the fundamental variables of the functional derivatives, rather than $\gamma = |\nabla\rho|^2$.

**Rationale**

The reason is given in Step 2 of Section 2.2 of [concept.md](concept.md): $\nabla\rho$ is a first-order quantity in the density matrix, while $\gamma$ is a second-order one. Using $\rho_r$ makes the subsequent program derivations considerably simpler.

The price is that the grid dimension increases from 1/2/3 (LDA/GGA/mGGA) to 1/4/5 (spin-unpolarized); but this increase happens at the eval_xc level and does not show up in the GEMM bottleneck.

In the implementation, LibXC uses $\gamma$ as its variable; the transformation from $\gamma$ derivatives to $\rho_r$ derivatives (sigma unfolding) is done in `xceff/xc_deriv.rs` via the chain rule:

$$
\frac{\partial(f\rho)}{\partial\rho_r} = 2 f^\gamma \rho_r, \quad r \in \{x, y, z\}
$$

For second- and third-order derivatives, the diagonal correction terms of $\partial^2/\partial\gamma^2$ must also be handled (in `transform_xc_inner`, adding the $2 f^{\gamma\gamma}$-related contributions to the diagonal elements with $\chi, \chi' \in \{x, y, z\}$).

## 6. Restrictions of LAPL-Type mGGA

**Decision**

- Density-on-grids generation supports the LAPL component (`XCDenType::LAPL`), but potential-matrix assembly and functional evaluation do not support LAPL.

**Rationale**

The current program can compute the LAPL density-on-grids component:

$$
\nabla^2 \rho_g = 4 \tau_g + 2 \sum_{\mu\nu} \varphi_{g \mu} D_{\mu\nu} (\varphi_{g \nu, xx} + \varphi_{g \nu, yy} + \varphi_{g \nu, zz})
$$

However, `contract_ao_wv` does not support LAPL contraction, and `libxc_eval_eff` raises an error directly for functionals that need the Laplacian. The reason is that LAPL contraction requires second-order AO derivatives ($n_\text{comp} = 10$), increasing the GEMM count, while LAPL-type functionals have limited application scenarios.

## 7. The Progressive Design of the XCDenType Enum

**Decision**

- The `XCDenType` enum (`RHO`/`SIGMA`/`TAU`/`LAPL`) is designed progressively: each higher-level type contains all components of the lower-level types.

**Rationale**

This allows `XCDenType` to control, simultaneously:

- the component count of the output density on grids, `num_nvar()` (1/4/5/6);
- the required AO derivative order, `num_ao_deriv()` (0/1/1/2);
- the required AO component count, `num_ao_comp()` (1/4/4/10, looked up in `AO_DERIV_DIM = [1, 4, 10, 20, 35]`).

The density-component order is unified as $\rho, \rho_x, \rho_y, \rho_z, \tau, \nabla^2\rho$, consistent with the progressive relation of `XCDenType`. Regardless of whether the functional needs $\tau$, $\tau$ is always the 5th component in the LAPL type (not the 4th), guaranteeing the consistency of component indices.

## 8. Grid Batching (nbatch) vs Grid Chunking (nchunk)

**Decision**

- `NIMatmul` uses `nbatch` to control memory batching and `nchunk` to control parallel chunking, set independently; `nbatch` should be a multiple of `nchunk`.

**Rationale**

`NIMatmul` has two granularity parameters:

- `nbatch`: the memory-control parameter. The full AO tensor `[ngrids, nao, ncomp]` can be too large (several GB for large systems), so the grid is processed in batches of `nbatch`. Each batch independently computes the AO, density, functional, and potential matrices, accumulating the results. The default is $1536 \times 1 \times n_\text{threads}$ (the thread count is determined by rayon at runtime).
- `nchunk`: the parallel-granularity parameter. In the pure-function layer, the grids are distributed to threads in chunks of `nchunk`. `nchunk` should match the KC dimension of GEMM (usually 256-512, for good cache utilization); the default 1536 is a few times KC.

The relation is full-grid > batch > chunk > per-grid = 1.

## 9. Grid-Space Contraction of the fxc/kxc Effective Potentials

**Decision**

- Functional outputs of second order and above do not enter potential assembly directly; they are first contracted in grid space with the perturbation densities into a first-order-form "effective potential", which then reuses the same `contract_ao_wv` contraction functions.

**Rationale**

For the second-order (fxc) and third-order (kxc) responses, the functional outputs are high-dimensional tensors ($f^{\chi\chi'}$ or $f^{\chi\chi'\chi''}$) and cannot be passed into `contract_ao_wv` directly. They must first be contracted in grid space:

- fxc: $\text{fxc\_eff\_contracted}^{\chi} = \sum_{\chi'} f^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}]$, yielding an effective potential of `[ngrids, nvar]`
- kxc: $\text{kxc\_eff\_contracted}^{\chi} = \sum_{\chi'\chi''} f^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'] \, \xi_{g}^{\chi''}[\mathbf{R}'']$, yielding an effective potential of `[ngrids, nvar]`

The contracted effective potential has the same shape as vxc_eff, so the `contract_ao_wv` function can be reused. This design keeps the interface of `contract_ao_wv` uniform across all derivative orders.

## 10. Spin-Dimension Conventions of UKS

**Decision**

- In open-shell tensors, the spin dimension immediately follows the density-variable dimension, and sits inside the grid dimension.

**Rationale**

The open-shell (UKS) tensor-shape conventions are:

- rho: `[ngrids, nvar, 2]` (not `[ngrids, 2, nvar]`)
- fxc_eff: `[ngrids, nvar, 2, nvar, 2]`
- fxc output: `[nao, nao, 2, nset]` (not `[nao, nao, nset, 2]`), which differs from PySCF.

The benefit of this convention is that the shapes of the functional outputs of all orders are constructed regularly (`nvar`, `2` repeated per order), and the component slicing by `XCDenType` (e.g. `rho.i((.., 1..4))`) is not disturbed by the spin dimension.

## 11. LDA Shape Convention

**Decision**

- The density on grids for LDA keeps its `nvar = 1` dimension, without squeezing.

**Rationale**

LDA has `nvar = 1`, i.e. only one density component. PySCF often squeezes this component away, making rho (in the closed shell) of shape `[ngrids]` rather than `[ngrids, 1]`.

But to keep the interface as uniform as possible, and to preserve the shape structure across all density types, we choose not to squeeze the LDA component, so that rho for LDA in the closed shell has shape `[ngrids, 1]`.

## 12. Optimization Paths and Deferred Items

**Decision**

- Performance optimization proceeds layer by layer by scope of impact: hotspots in the pure-function layer, data structures in the public-interface layer, and a stable driver layer; some optimizations are explicitly deferred.

**Rationale**

The core goal of this module is API-design concepts and a correct implementation, not a performance showcase; but the separation of the three layers provides a clear path for subsequent optimization. The discussion below is classified by the scope of changes each optimization requires.

Hotspots directly optimizable at the bottom layer:

- **Grid sparsity (non0tab)**: the current program stores dense $\phi_{g \mu}$. Introducing PySCF-style non0tab sparse masks would require adding a mask field to `NIMatmul`, generating the mask in `prepare_ao`, and adding mask parameters or new pure functions in the pure-function layer. This affects the middle layer and the pure-function layer, but not the driver layer.
- **Psi4-style blocking**: compressing $\phi_{g \mu}$ into $\phi_{g \mu'}^{\text{packed}}$ plus a mapping table would require new data structures in the middle layer and corresponding contraction functions in the pure-function layer. Likewise no impact on the driver layer.
- **The `contract_ao_wv` family** is the computational bottleneck. The current implementation uses standard GEMM. Considering the sparsity of DFT grids and basis functions, and that the grid count usually far exceeds the AO count, dedicated micro-kernels for grid-basis products could be developed in the future, to better exploit cache and SIMD instructions. For small systems GEMM is already very efficient; for large systems, customized kernels may bring significant gains.

Optimizations requiring coordination across layers:

- **Bra-ket low-rank optimization**: `homogeneous_braket` and `one_bra_mult_ket` forms are already supported. For the bra-trans variant of fxc (`make_fxc_pot_with_eff_bra_trans`), the output changes from `[nao, nao]` to `[nao, nocc]`, which affects the data structures of the upper-level algorithms (such as TDDFT solvers). Therefore, generalizing bra-trans requires joint adaptation of the driver layer and its callers.
- **Grid batching strategy**: the choice of `nbatch` affects memory usage and parallel efficiency. For particularly large systems, finer batching strategies may be needed (e.g. batching by atom rather than by grid order), which requires adjusting the batching logic of `NIMatmul` in the middle layer.
- **Functional-evaluation parallelism**: the default chunk sizes of `libxc_eval_eff_parallel` are adjusted by density type. Under nested parallelism (the driver layer already inside a rayon thread pool), it must degrade to serial. More flexible parallel strategies (e.g. an independent thread pool) may be needed in the future.

Optimizations deferred for now:

- **Exploiting $\mu, \nu$ symmetry (SYRK-type optimizations)**: the computations of $\rho_g$ and $\tau_g$ can exploit the density-matrix symmetry to halve the cost, but $\rho_g^t$ cannot. Considering that GGA/mGGA are the mainstream workloads, and that SYRK micro-kernels differ from GEMM ones, this is not a priority for now (see Section 3).
- **Complex-number support**: currently only `f64` is supported. Complex types would require new micro-kernels and conjugation handling; not implemented for now.
- **LAPL-type mGGA potential contraction**: see Section 6.
