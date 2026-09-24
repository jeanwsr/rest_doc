# API Design Concepts for DFT Grid Integration

This document discusses several core concepts in the API design of DFT grid integration. We hope these concepts will help the future design of DFT grid integration in terms of extensibility, performance optimization, and API friendliness.

This document also attempts to clarify the terms, concepts, and formulas that a DFT program needs, in a relatively consistent text. The statements here are not necessarily the most orthodox, but they are indeed the line of thought consistently followed when designing the API.

:::{note}
The formula notation and tensor shapes of this document follow the conventions of [def.md](def.md): $\chi$ denotes the component index of the basic DFT variables $\boldsymbol{\xi}$; the spatial index is denoted by $r$ in density expressions and by $t$ in response contractions. Tensor shapes follow the column-major convention.

This document only treats the closed-shell (RKS) formulas; the shape differences of the open-shell case are given in Section 4 of def.md and in [basic-api.md](basic-api.md). The function interfaces involved in this document are listed in [basic-api.md](basic-api.md), and the design trade-offs and their reasons are given in [adr.md](adr.md).
:::

## 1. Generalities: Key API Design for Self-Consistent Field

### 1.1 Premise: Energy Decomposability and Its Relation to the Density Matrix

The fundamental variables of the computational task are basis functions; these basis functions can be atomic orbitals, plane waves, real-space grids, though atomic orbitals are the most common. In this project, basis functions are customarily called AO (atomic orbital), but this does not mean we can only handle atomic orbitals.

**The energy is decomposable**. For example, for the wB97X-V functional, its energy can be split into the following parts:

$$
E[\mathbf{D}] = E_\text{nuc-repl} + E_\text{kin} + E_\text{nuc} + E_\text{J} + E_\text{K} + E_\text{srK} + E_\text{xc} + E_\text{VV10}
$$

One key point is that, assuming nuclei are not treated as variables, **the energy and its components are functions of the density matrix $\mathbf{D}$**. This statement itself is simple, but it implies the following corollaries or less-common counterexamples:

- The density matrix is defined as a function of the basis-function coefficient matrix $C_{\mu i}$:

    $$
    D_{\mu\nu} = \sum_i C_{\mu i} C_{\nu i}
    $$

    Therefore, the energy and its components can also be viewed as functions of the coefficient matrix $\mathbf{C}$.

- The energy cannot be a function of the coefficient matrix $\mathbf{C}$ alone. It must be expressible as an explicit function of the density matrix $\mathbf{D}$. The counterexamples here are the various orbital-optimized post-SCF methods.

- The relation between the energy and the density must strictly satisfy the Hartree-Fock-Roothaan equations:

    $$
    \mathbf{V}[\mathbf{D}] \mathbf{C} = \mathbf{S} \mathbf{C} \boldsymbol{\epsilon}
    $$

    where the Fock matrix $\mathbf{V}$ is the gradient of the energy with respect to the density matrix $\mathbf{D}$:

    $$
    \mathbf{V}[\mathbf{D}] = \frac{\partial E[\mathbf{D}]}{\partial \mathbf{D}}
    $$

    The counterexamples here are methods such as constraint DFT (cDFT) or density-corrected DFT (dcDFT). cDFT contains constraint terms; the non-self-consistent character of dcDFT makes its Fock matrix no longer the gradient of the energy. Although we can still make use of the current API when writing programs, in terms of program design philosophy, and to keep the discussion simple, we do not consider these methods for now.

- Atomic structure is out of scope; we only consider electronic structure (the electron density represented in basis functions). Therefore, for the nuclear repulsion energy $E_\text{nuc-repl}$, as well as methods similar to DFT-D3, we regard their energy terms as constants.

Also note that **the Fock matrix in this document is denoted $\mathbf{V} [\mathbf{D}]$, not the usual $\mathbf{F} [\mathbf{D}]$**. The next subsection makes the notation of this document explicit.

### 1.2 The Core Problem: Derivatives of the Energy with Respect to the Density Matrix

Many of the core technical problems in calculations involving self-consistent-field methods boil down to how to compute derivatives of the energy with respect to the density matrix.

We need the following definitions:

| Derivative order | Variable | Full notation | Customary program symbols |
|--|--|--|--|
| 1 | $\mathbf{V}$ | $\mathbf{V} [\mathbf{D}]$ | `vxc`, `v`, `fock`, `veff` |
| 2 | $\mathbf{F}$ | $\mathbf{F} [\mathbf{D}, \mathbf{R}]$ | `fxc`, `f`, `resp` |
| 3 | $\mathbf{K}$ | $\mathbf{K} [\mathbf{D}, \mathbf{R}^1, \mathbf{R}^2]$ | `kxc`, `k` |

- In the nimatmul program, the three quantities above are named `vxc`, `fxc`, `kxc` according to the functional-derivative order, corresponding to the XC grid-integration contributions to $\mathbf{V}$, $\mathbf{F}$, $\mathbf{K}$ respectively.

- First-order derivative definition

    $$
    V_{\mu \nu} = \frac{\partial E}{\partial D_{\mu \nu}}
    $$

    The first-order derivative is usually denoted $F_{\mu \nu}$ in other literature. Here we deliberately use a different symbol to distinguish first- and second-order derivatives.

- Second-order derivative definition

    $$
    F_{\mu \nu} = \sum_{\kappa \lambda} \frac{\partial^2 E}{\partial D_{\mu \nu} \, \partial R_{\kappa \lambda}} R_{\kappa \lambda}
    $$

    where $R_{\kappa \lambda}$ is a perturbation density matrix; in TDDFT (Casida equations) it is often an excited-state density matrix, and in gradient problems (CP-HF/KS equations) it is a derivative density matrix.

    In other literature, the second-order derivative is usually written $\sum_{\kappa \lambda} A_{\mu \nu, \kappa \lambda} R_{\kappa \lambda}$, or $G_{\mu \nu} [\mathbf{R}]$.

- Third-order derivative definition

    $$
    K_{\mu \nu} = \sum_{\kappa \lambda} \sum_{\kappa' \lambda'} \frac{\partial^3 E}{\partial D_{\mu \nu} \, \partial R^1_{\kappa \lambda} \, \partial R^2_{\kappa' \lambda'}} R^1_{\kappa \lambda} R^2_{\kappa' \lambda'}
    $$

    where $R^1_{\kappa \lambda}$ and $R^2_{\kappa' \lambda'}$ are two perturbation density matrices.

All of the definitions above apply to every energy component usable in self-consistent-field calculations, and they can be linearly summed. For instance, for the wB97X-V functional, its Fock matrix can be written, in the same manner, as

$$
\mathbf{V} [\mathbf{D}] = \mathbf{V}_\text{nuc-repl} + \mathbf{V}_\text{kin} + \mathbf{V}_\text{nuc} + \mathbf{V}_\text{J} + \mathbf{V}_\text{K} + \mathbf{V}_\text{srK} + \mathbf{V}_\text{xc} + \mathbf{V}_\text{VV10}
$$

Therefore, for any energy component, one can define the following program interface (shown as Python pseudocode):

```python
class EnergyComponent:
    def energy(self, dm0: np.ndarray) -> float:
    def get_fock(self, dm0: np.ndarray) -> np.ndarray:
    def get_2nd_resp(self, dm0: np.ndarray, dm1: np.ndarray) -> np.ndarray:
    def get_3rd_resp(self, dm0: np.ndarray, dm1: np.ndarray, dm2: np.ndarray) -> np.ndarray:
```

### 1.3 Technical Considerations: Performance and Extensibility

We believe the interface design above can realize all important computational-chemistry needs. In practice, however, we also need the following considerations and adaptations:

- **Closed-shell and open-shell density matrices differ in shape**. The closed-shell shape is `[nao, nao]`, while the open-shell shape is `[nao, nao, 2]` (col-major).

- Density matrices, especially perturbation density matrices, **often come in multiples** rather than singly. Therefore, for the closed shell, the passed `dm1` and `dm2` should generally be allowed to be 3-dim tensors, or lists of 2-dim matrices. The passed `dm0` is in most cases a single matrix.

- Density matrices are usually low-rank. They are built from occupied-orbital coefficients, and the occupation $n_\text{occ}$ is usually far smaller than the basis size $n_\text{AO}$. Even the perturbation densities required for solving the CP-KS equations (which differ from the self-consistent-field density) can, by construction (or via decompositions such as SVD), be factored into products of two $n_\text{AO} \times n_\text{occ}$ matrices. Therefore **the low-rank structure of density matrices should be exploited as much as possible** to reduce the computational cost.
    - On the other hand, for particularly large molecules, the density matrix itself may be sufficiently sparse. The sparsity of zero values in a density matrix and its low-rank structure cannot (or at least can hardly) be exploited simultaneously; in some cases the sparsity of the density matrix itself can also be exploited. But for atomic systems of no more than 100 atoms, the low-rank structure of the density matrix is usually more important.
    - From the API design perspective, exploiting the low-rank property is equivalent to passing occupied-orbital coefficient matrices into the computational functions:

        ```python
        def get_fock_by_occ(self, occ_coeff: np.ndarray) -> np.ndarray:
        ```

        The concrete interface form may vary somewhat, but the key idea is that **the API design should allow passing occupied-orbital coefficient matrices directly**, to exploit the low-rank structure of density matrices. We will see later how this is concretely exploited in DFT grid integration.

## 2. Computational Steps of DFT Grid Integration

### 2.1 Comparison with Other Energy Components

DFT grid integration shares similarities with, and has differences from, other energy-component contributions in computational chemistry.

The similarities are

- All the general concepts above apply to DFT grid integration.
- DFT likewise involves the question of whether occupied orbitals can be exploited for optimization. Introducing occupied orbitals affects not only implementation details, but also complicates the API design problem.
- The computational complexity of DFT grid integration is relatively low, but the workload is often not small, especially for small to medium systems. Therefore, performance still needs to be considered when designing the API.

The differences are

- For J/K computations, on the one hand the energy is a quadratic function of the density matrix (so the J/K energy has no third-order term $\mathbf{K}$), and on the other hand the second-order J/K computation closely resembles the first-order one. Therefore, although J/K is computationally heavy and approximating it remains one of the research focuses of electronic-structure programs, its API design is comparatively simple.
- The difficulty of DFT lies in that its energy is a function of the density (it cannot be Taylor-truncated to finite order). Moreover, the higher the order, the larger the workload. Without reasonable program design, its implementation can become very hard.

### 2.2 Density Generation, Functional Evaluation, and Assembly of Energy-Derivative Matrices

The computation of DFT grid integration can be clearly split into three steps. These steps have clear input-output dependencies, but differ significantly in computational character. We describe them one by one for the closed-shell (RKS) case, distinguishing the three density types $\text{RHO}$ (LDA), $\text{SIGMA}$ (GGA), and $\text{TAU}$ (mGGA) (as in def.md, we do not consider $\text{LAPL}$-type mGGA for now; the current status of LAPL is given in [adr.md](adr.md)).

#### Step 1: Density-on-grids generation (eval_rho)

Taking the density matrix $D_{\mu\nu}$ as input, with the AO on grids $\phi_{g \mu}$ and their spatial derivatives $\phi_{g \mu}^r$ (where $r \in \{x, y, z\}$) as intermediates, generate the density on grids. Define the density variables

$$
\xi_g[\mathbf{D}] := (\rho_g, \rho_g^x, \rho_g^y, \rho_g^z, \tau_g)
$$

For the closed shell, this variable is denoted $\xi_g^{\chi}$ and is a 2-dim array `[ngrids, nvar]` in the program. For the open shell, this variable is denoted $\xi_g^{\chi \sigma}$, i.e. with an additional spin dimension, and is a 3-dim array `[ngrids, nvar, 2]` in the program.

The component formulas are as follows:

- **Density** $\rho_g$ (needed by all types)

    $$\rho_g = \sum_{\mu\nu} \phi_{g \mu} D_{\mu\nu} \phi_{g \nu}$$

    The algorithm is one matrix multiplication plus one scaled reduction:

    $$\bar{\phi}_{g \mu} = \sum_\nu D_{\mu\nu} \phi_{g \nu} \quad \text{(GEMM)}$$

    $$\rho_g = \sum_\mu \phi_{g \mu} \bar{\phi}_{g \mu} \quad \text{(memory bounded)}$$

- **Density gradient** $\rho_g^r$ (needed by GGA / mGGA)

    $$\rho_g^r = 2 \sum_{\mu\nu} \phi_{g \mu}^r D_{\mu\nu} \phi_{g \nu}$$

    Since $\bar{\phi}_{g \mu}$ is already obtained in the computation of $\rho_g$, each gradient component needs only one additional scaled reduction $\sum_\mu \phi_{g \mu}^r \bar{\phi}_{g \mu}$, multiplied by a factor of 2. The extra cost of the three components together is memory bounded.

- **Kinetic energy density** $\tau_g$ (needed by mGGA)

    $$\tau_g = \sum_{r,\mu\nu} \frac{1}{2} \phi_{g \mu}^r D_{\mu\nu} \phi_{g \nu}^r$$

    This requires one new matrix multiplication per $r$, $\bar{\phi}_{g \mu}^{(r)} = \sum_\nu D_{\mu\nu} \phi_{g \nu}^r$ (with factor 1/2), followed by the scaled reduction $\sum_\mu \phi_{g \mu}^r \bar{\phi}_{g \mu}^{(r)}$. The three GEMMs cost about $6 n_\text{basis}^2 n_\text{grid}$ FLOPs.

**Low-rank optimization**: when the density matrix has the low-rank structure $D_{\mu\nu} = \sum_i C_{\mu i} C_{\nu i}$ ($n_\text{occ} \ll n_\text{basis}$), the cost can be reduced through the occupied-orbital coefficients $C_{\mu i}$ (bra-ket form):

$$\phi_{g i} = \sum_\mu \phi_{g \mu} C_{\mu i} \quad \text{(GEMM)}$$

$$\rho_g = \sum_i \phi_{g i}^{(\text{bra})} \phi_{g i}^{(\text{ket})} \quad \text{(memory bounded)}$$

For 3-$\zeta$ basis sets of the def2-TZVP class, $n_\text{basis} / n_\text{occ} \sim 10$, so the cost drops by roughly an order of magnitude.

The cost of density-on-grids generation for each density type is summarized below:

| Density type | Variable count $n_\text{var}$ | AO derivative order | AO component count $n_\text{comp}$ | GEMM count (DM input) |
|--|--|--|--|--|
| RHO (LDA) | 1 | 0 | 1 | 1 |
| SIGMA (GGA) | 4 | 1 | 4 | 1 |
| TAU (mGGA) | 5 | 1 | 4 | 4 |

Density-on-grids generation is the step second only to response-matrix generation in cost.

#### Step 2: Functional evaluation (eval_xc)

Substitute the density on grids $\xi_g[\mathbf{D}]$ into the density functional $f(\rho, \nabla\rho, \tau)$ and its partial derivatives of each order, obtaining the functional outputs. This part costs $O(n_\text{grid})$ and is not the bottleneck.

The output dimensions of the functional are defined as follows (closed shell):

- First-order derivative (effective potential, for vxc): $f_{g}^{\chi}$, shape `[ngrids, nvar]`
- Second-order derivative (effective kernel, for fxc): $f_{g}^{\chi \chi'}$, shape `[ngrids, nvar, nvar]`
- Third-order derivative (effective kernel, for kxc): $f_{g}^{\chi \chi' \chi''}$, shape `[ngrids, nvar, nvar, nvar]`

For the open shell, a spin dimension of 2 is inserted after each `nvar`; for example the first-order derivative is `[ngrids, nvar, 2]`.

:::{note}
**Important design choice**: this program uses the density-gradient components $\rho_g^r$ as the fundamental variables of the functional, rather than $\gamma = |\nabla\rho|^2$. The reason is that $\gamma$ is a second-order quantity in the density matrix (taking an additional density-matrix derivative through $\gamma$ does not vanish), whereas $\nabla\rho$ is strictly a first-order quantity in the density matrix. This makes the subsequent formula derivations and program implementation considerably simpler. The price is larger grid dimensions (for spin-unpolarized LDA/GGA/mGGA the variable count increases from 1/2/3 to 1/4/5); since the bottleneck of DFT grid integration is GEMM rather than grid dimensions, this price usually does not show up in the real bottleneck. A more detailed discussion is given in [adr.md](adr.md).
:::

The transformation from $\gamma$ to $\rho_r$ follows the chain rule:

$$
\frac{\partial(f\rho)}{\partial\rho_r} = \frac{\partial(f\rho)}{\partial\gamma} \frac{\partial\gamma}{\partial\rho_r} = 2 f^\gamma \rho_r, \quad r \in \{x, y, z\}
$$

#### Step 3: Assembly of energy-derivative matrices (contract_ao_wv)

Contract the grid-wise weighted quantity $w_g f_\text{eff}$ with the AO on grids $\phi_{g \mu}$, obtaining the energy-derivative matrix (the XC contribution to the Fock matrix).

**First-order response (vxc)**. For the closed shell:

$$V_{\mu\nu}^{\text{xc}} = \sum_g w_g f_g^\rho \, \phi_{g \mu} \phi_{g \nu} \quad \text{(LDA)}$$

$$V_{\mu\nu}^{\text{xc}} \leftarrow \sum_{t,g} w_g f_g^{\rho_t} \phi_{g \mu}^t \phi_{g \nu} + \text{swap}(\mu, \nu) \quad \text{(GGA)}$$

$$V_{\mu\nu}^{\text{xc}} \leftarrow \sum_{t,g} \frac{1}{2} w_g f_g^\tau \phi_{g \mu}^t \phi_{g \nu}^t \quad \text{(mGGA)}$$

The common structure of the expressions above is: **bra $\phi_{g \mu}^{(\text{lhs})}$, ket $\phi_{g \nu}^{(\text{rhs})}$, and grid-wise weighted quantity $w_g f_\text{eff}$, finally summed over the grid index $g$**. The bra, ket, and weighted quantity differ in content across density types, but the structure is the same. In the implementation, the asymmetric half matrix is computed first and symmetrized in one step; the contraction coefficients differ across density types, see [adr.md](adr.md).

**Second-order response (fxc)**. Introducing the basic DFT variables $\boldsymbol{\xi}$, the fxc expression is

$$
F_{\mu\nu}[\mathbf{R}] = \sum_g \sum_\chi w_g \left(\sum_{\chi'} f_g^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}]\right) \frac{\partial \xi_g^{\chi}[\mathbf{D}]}{\partial D_{\mu\nu}}
$$

where $\sum_{\chi'} f_g^{\chi\chi'} \xi_{g}^{\chi'}[\mathbf{R}]$ is a contraction in grid space, yielding a grid quantity indexed by $\chi$; the subsequent $\sum_\chi (\cdots) \, \partial \xi_g^{\chi} / \partial D_{\mu\nu}$ is handled exactly as in vxc. This means the fxc grid-contraction step can reuse the vxc `contract_ao_wv` function.

**Third-order response (kxc)** similarly, after contracting in grid space $\sum_{\chi'\chi''} f_g^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'] \, \xi_{g}^{\chi''}[\mathbf{R}'']$, `contract_ao_wv` is reused as well.

**Computational bottleneck**. Response-matrix generation is the most costly step of DFT grid integration. For vxc, the FLOPs scale as $O(n_\text{basis}^2 n_\text{grid} n_\text{var})$; for fxc, with $n_\text{set}$ perturbation density matrices, it scales as $O(n_\text{basis}^2 n_\text{grid} n_\text{var}^2 n_\text{set})$. Density-on-grids generation comes second. Functional evaluation is negligible.

### 2.3 Overall Design: Separating Functional Evaluation from Grid Contraction

From the analysis of the previous section, the three computational steps differ significantly in character:

- **Density-on-grids generation** and **assembly of energy-derivative matrices** are both matrix operations (GEMM + scaled reduction), independent of the specific functional;
- **Functional evaluation** is a per-grid-point operation ($O(n_\text{grid})$), independent of the AO basis structure.

Therefore, we separate functional evaluation and grid contraction (called NIMatmul in this program) into independent modules. This separation brings the following benefits:

1. **The grid-contraction functions can take "effective potentials" (eff_pot) as input**, rather than the raw functional outputs. An effective potential is the grid-space contraction of the functional output $f_g^{\chi}$ with the density on grids $\xi_{g}^{\chi'}[\mathbf{R}]$ (for second order and above), yielding a grid-wise weighted vector. This makes the `contract_ao_wv` family of functions completely independent of LibXC, testable and optimizable in isolation.

2. **The functional-evaluation module can be replaced independently**. LibXC is used at present, but in the future XCFun, machine-learned functionals, or other functional engines can be plugged in, as long as their output format conforms to the effective-potential convention.

3. **Clear data flow**. The complete data flow is:

    ```
    // vxc
    dm → [eval_rho] → rho → [eval_xc] → vxc_eff → [contract_ao_wv] → vxc

    // fxc
    dm + dm1 → [eval_rho] → rho, rho1 → [eval_xc] → fxc_eff → [contract with rho1] → fxc_eff_contracted → [contract_ao_wv] → fxc

    // kxc
    dm + dm1 + dm2 → [eval_rho] → rho, rho1, rho2 → [eval_xc] → kxc_eff → [contract with rho1, rho2] → kxc_eff_contracted → [contract_ao_wv] → kxc
    ```

    For the first order (vxc), the functional output is directly the effective potential. For the second order (fxc) and third order (kxc), one must first contract in grid space $\sum_{\chi'} f^{\chi\chi'} \xi_{g}^{\chi'}[\mathbf{R}]$ (or $\sum_{\chi'\chi''} f^{\chi\chi'\chi''} \xi_{g}^{\chi'}[\mathbf{R}'] \xi_{g}^{\chi''}[\mathbf{R}'']$), and only then enter `contract_ao_wv` with the effective potential.

4. **Low-rank optimization can be implemented at the eval_rho level**, without affecting the interfaces of eval_xc and contract_ao_wv. This program provides four ways of density-on-grids generation (see [basic-api.md](basic-api.md)), supporting low-rank optimization in different scenarios.

### 2.4 Three-Layer Architecture: Pure-Function Algorithms, DFT Public Interface, Driver Interface

This program adopts a three-layer architecture. From bottom to top:

#### Pure-function layer (bottom)

This layer corresponds, in REST, to the functions in `numint_matmul/pure_eval_rho.rs` and `numint_matmul/pure_xcpot.rs`. Their characteristics:

- **Stateless**: functions hold no `self` or hidden state; all inputs are passed explicitly as parameters.
- **Simple parameters**: inputs are tensor views and flag enums (e.g. `XCDenType`); outputs are written into pre-allocated buffers (`*_with_output` suffix). Complicated types (e.g. a complete grids structure) are best avoided. Ideally, such functions are also easy to export to a C API.
- **Flat data structures**: no complicated logic such as grid batching or AO caching.

The pure functions are also where the **performance hotspots** live. Since the parameter lists are fully explicit, this layer can be tested, replaced, or further optimized independently, without affecting the upper interfaces.

#### DFT public-interface layer (middle)

This layer is the `NIMatmul` struct (`numint_matmul/nimatmul.rs`) and the functional-evaluation bridge (the `dft/xceff` module: `flags.rs`, `libxc_wrap.rs`, `xc_deriv.rs`). Their characteristics:

- **Stateful**: `NIMatmul` manages the AO cache, grid batching parameters (`nchunk`, `nbatch`), the integral engine, and so on.
- **More complex parameters**: method calls involve cache strategies, batching logic, functional-engine selection, and so on.

This is the layer **where changes are likely to happen**. For example, replacing the integral engine, changing the grid storage format (from dense $\phi_{g \mu}$ to Psi4 blocking or sparse formats), or adjusting cache strategies are all mainly implemented in this layer. Importantly, however, **changes in this layer should not affect the pure-function layer and the driver layer**.

#### Driver layer (top)

This layer chains the computational steps into complete workflows, facing the end user. In REST, this layer is the Hessian driver interface (`numint_matmul/hess_rks.rs`, `hess_uks.rs`, including the `RHessKSNIMatmul` / `UHessKSNIMatmul` structs and the traits they implement); its role is analogous to the driver functions of PySCF's `dft.numint.nr_rks` / `nr_uks`. As long as the pure-function layer and the public-interface layer keep their APIs, the driver layer can remain stable.

**The core meaning of the three-layer architecture** is that performance optimization should concentrate on the pure-function layer (hotspot functions such as `contract_ao_wv`), data-structure adjustments and non-trivial performance problems happen in the middle layer, while the top-level user interface stays unchanged. This allows development work at different levels to proceed in a decoupled manner.

In the future, we also plan to abstract the DFT layer and the driver layer into traits, further allowing various data structures to implement the same interface and enhancing flexibility and extensibility. In REST's Hessian driver layer, this idea has been partially realized: `RHessKSNIMatmul` / `UHessKSNIMatmul` plug into the Hessian solver through traits such as `HessUtilAPI` and `RHessElecInteractAPI` (see [skeleton2.md](skeleton2.md) and [vmat1.md](vmat1.md)).
