# TDDFT AO mode: transition-density kernel

This page describes the matrix-vector products in `tddft_mode = "ao"`. The notation follows the [RI documentation](../ri.md) and [tddft-mo](tddft-mo.md): $\mu\nu$ atomic orbitals, $i,j$ occupied MOs, $a,b$ virtual MOs, $P$ auxiliary basis, $g$ grid point, $\alpha,\beta$ density-variable components, $\mathbb{A}$ the trial-vector set index; underlined subscripts (e.g. $\underline{g}$) denote batched indices.

## Overview

The MO mode must pre-store three MO-basis RI tensors totaling $n_\mathrm{aux}(n_\mathrm{occ}n_\mathrm{vir} + n_\mathrm{occ}^2 + n_\mathrm{vir}^2)$; when $n_\mathrm{occ}, n_\mathrm{vir} \gg n_\mathrm{basis}$ this far exceeds the AO-basis $n_\mathrm{aux} n_\mathrm{basis}^2$ (i.e. `scf.rimatr`), and the MO transformation itself is a one-off large DGEMM. The AO-mode strategy is:

- the Davidson iteration stays in the MO amplitude space (small subspace);
- each matrix-vector product transforms the whole trial-vector block into **AO transition densities** $P_{\mu\nu}^{\mathbb{A}}$;
- Coulomb and exchange directly call the density-driven interface of `ri_jk` (reusing the SCF's `scf.rimatr`; no new large tensors);
- the XC kernel is evaluated on the numerical grid for the whole trial-vector block in one batched pass via `dft::numint_matmul::NIMatmul` (following the `vind(zs)` design of PySCF's `_gen_tda_operation`);
- the result is contracted back into the amplitude space with the MO coefficients.

## Function `prepare_ao_data_with_spin`

Function path: `ri_tddft::tddft::prepare_ao_data_with_spin` (signature `(scf, tddft_spin: Option<&str>)`; `tddft_solver` passes the deck-resolved spin channel, while the stability module passes `"singlet"`/`"triplet"` directly, ignoring the deck)

Builds the AO-mode `TDDFTData`. The MO-specific members (`fxc`, `fxc_u`, `ri_terms`) are empty/`None`; the AO-side members are (one entry per spin sector for an unrestricted reference, `Vec` length = `n_sectors()`):

| Member (`TDDFTData`) | Meaning | Dimensions | Remarks |
|--|--|--|--|
| `c_occ` | $C_{\mu i}$ | $(n_\mathrm{basis}, n_\mathrm{occ})$ | per sector; empty sectors (occ = 0) carry zero-column matrices |
| `c_vir` | $C_{\mu a}$ | $(n_\mathrm{basis}, n_\mathrm{vir})$ | per sector |
| `ni` | `NIMatmul` numerical integrator | — | libcint AO cache; real grid weights |
| `fxc_eff` | raw kernel table (singlet factor included) | restricted $(n_\mathrm{grid}, n_\mathrm{var}, n_\mathrm{var})$; unrestricted $(n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2)$ | weights not multiplied; see below |
| `den_type` | `RHO` / `SIGMA` | — | determines $n_\mathrm{var}$ |
| `grid_batch` | grid-batching switch | — | |
| `fxc_driver` | `Option<FxcDriver>` | — | one of `"mo"`/`"semitrans"`/`"dm"`; `None` only for an HF reference (no kernel, J/K only) |
| `psi_occ` | occ-MO projections on grid $\psi_i(g)$ | $[n_\mathrm{grid}, n_\mathrm{occ}]$ | built only for the `MO`/`SEMITRANS` fxc drivers; projected grid-batch-wise, t-ready layout |
| `psi_occ_grad` | $\partial_d \psi_i(g)$ | $[3, n_\mathrm{grid}, n_\mathrm{occ}]$ | GGA only; leading $d$ axis keeps batch slices contiguous |

Two early-exit branches:

- **HF reference** (no libxc components): `fxc_eff`/`ni`/`fxc_driver` are all `None`, and the matrix-vector products run the RI J/K parts only. This lets AO-mode TDDFT and the stability analysis run for HF references without any DFT grid;
- **RSH functionals**: the response exchange requires the short-range three-center integrals `scf.rimatr_sr` (built by the SCF with the same functional); a missing table panics with a hint to re-run the SCF.

**Kernel-table convention**: `fxc_eff` stores the **raw (unweighted) kernel** with the singlet factor ($\times 2$) already applied, while `NIMatmul` is constructed with the real grid weights and multiplies them internally inside `make_fxc_pot_with_eff`. This differs from the MO-mode convention (weights pre-multiplied into `wfxc`), but the final contractions are mathematically equivalent.

**Spin-channel selection** (restricted references, driven by `tddft_spin`; the explicit argument of `prepare_ao_data_with_spin` overrides the deck):

| Channel | Kernel | Evaluation |
|--|--|--|
| Singlet | $f_s = 2 f_u$ | unpolarized |
| Unpolarized 'R' | $f_u$ | unpolarized (factor 1) |
| Triplet | $f_t = f_{\uparrow\uparrow} - f_{\uparrow\downarrow}$ | spin-polarized (see below) |
| `tddft_spin = "both"` | one kernel each for singlet and triplet | `tddft_main`'s `run_spin` calls this function per channel; the kernels are not shared |
| Unrestricted reference | spin-resolved $f_{\sigma_1\sigma_2}$ | spin-polarized evaluation at the real $(\rho_{\alpha 0}, \rho_{\beta 0})$; `fxc_eff: [n_\mathrm{grid}, n_\mathrm{var}, 2, n_\mathrm{var}, 2]`, no singlet/triplet factors |

The triplet kernel cannot be obtained from an unpolarized evaluation: `prepare_ao_data_with_spin` evaluates with `LibXCSpin::Polarized` at $\rho_\uparrow = \rho_\downarrow = \rho/2$ (halved gradients for GGA), giving the spin-resolved kernel $K[g, y_1, s_1, y_2, s_2]$ (REST's polarized transform already chain-rules to per-spin gradient components), then combines along the antisymmetric direction:

$$
f_t[y_1, y_2](g) = \frac{1}{2} \sum_{s_1 s_2} (\pm 1)^{s_1 + s_2}\, K[g, y_1, s_1, y_2, s_2]
$$

This is the triplet recipe of PySCF's `nr_rks_fxc_st` (CPL, 256, 454).

## Function `transition_density`

Function path: `ri_tddft::matvec_ao::transition_density`

$$
P_{\mu\nu}^{\mathbb{A}} = \sum_{ia} C_{\mu i}\, z_{ia}^{\mathbb{A}}\, C_{\nu a}
$$

One DGEMM. **Note that $P^{\mathbb{A}}$ is generally not symmetric** (the occupied and virtual index sets differ) — this is the central fact all subsequent AO-mode operators must handle; the XC kernel and the Coulomb term depend only on its symmetric part, while the exchange term must retain the full non-symmetric density.

## Function `get_j_ao_batched`

Function path: `ri_tddft::matvec_ao::get_j_ao_batched` (internally `ri_jk::pure_incore::get_vj_ri_incore_nonsym`)

Since $(\mu\nu|\kappa\lambda)$ is symmetric under $\kappa \leftrightarrow \lambda$, the Coulomb term depends only on the symmetric part of the transition density. The density is folded and routed through the standard packed contraction of the RI-JK incore algorithm (notation as in the `get_vj_ri_incore` section of the [ri documentation](../ri.md)):

$$
\begin{aligned}
D_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} &=
\begin{cases}
P_{\mu\nu}^{\mathbb{A}} + P_{\nu\mu}^{\mathbb{A}}, & \mu \neq \nu \\
P_{\mu\mu}^{\mathbb{A}}, & \mu = \nu
\end{cases}
&& \text{(eq.1 fold } \tilde\bowtie\text{)} \\
\mathscr{T}_{P}^{\mathbb{A}} &= \sum_{\mathrm{tp}(\mu\nu)} D_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} Y_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.2)} \\
J_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}} &= \sum_{P} \mathscr{T}_{P}^{\mathbb{A}} Y_{\mathrm{tp}(\mu\nu), P}
&& \text{(eq.3)} \\
J_{\mu\nu}^{\mathbb{A}} &\bowtie J_{\mathrm{tp}(\mu\nu)}^{\mathbb{A}}
&& \text{(eq.4)}
\end{aligned}
$$

The whole trial-vector block is packed as $[\mu, \nu, \mathbb{A}]$ in a single call. Unlike the SCF's `get_vj_ri_incore` (symmetric densities only, folds $2D - \mathrm{diag}$), the fold in `get_vj_ri_incore_nonsym` follows eq.1 above and is valid for non-symmetric transition densities.

## Function `get_k_ao_batched`

Function path: `ri_tddft::matvec_ao::get_k_ao_batched` (internally `ri_jk::pure_incore::get_vk_ri_incore_dm`)

The exchange term uses the density-driven incore algorithm of `ri_jk`, which **naturally handles non-symmetric densities** (per auxiliary column $M_{\mu\nu,\underline{P}} \cdot D \cdot M_{\nu\mu,\underline{P}}$, no folding step; equations in the `get_vk_ri_incore_dm` section of the [ri documentation](../ri.md)):

$$
K_{\mu\nu}^{\mathbb{A}} = \sum_{\kappa \underline{P}} \mathscr{T}_{\mu\kappa, \underline{P}}^{\mathbb{A}}\, Y_{\nu\kappa, \underline{P}} \qquad
\mathscr{T}_{\mu\kappa, \underline{P}}^{\mathbb{A}} = \sum_{\lambda} Y_{\mu\lambda, \underline{P}} D_{\kappa\lambda}^{\mathbb{A}}
$$

The **B block** has the exchange index ordering $(ib|aj)$, realized via the transposed density: $K\!\left[D^{\mathbb{A}^{\mathrm{T}}}\right]$, i.e. swapping $c_\mathrm{occ}/c_\mathrm{vir}$ and passing the transposed amplitudes.

**Driver selection** `tddft_ao_rik_driver` (AO mode only):

| Value | Path | Remarks |
|--|--|--|
| `"semitrans"` (default) | `get_vk_ri_incore_coeff_pair` | occupied-side semi-transformation then contract (below); **exact** at $O(n_\mathrm{occ})$ cost |
| `"dm"` | `get_vk_ri_incore_dm` | exact, whole-block batched, per-aux-column density-driven |
| `"lowrank"` | `get_vk_ri_incore_dm_lowrank` | per-vector SVD low-rank (threshold `tddft_svd_tol`), lossy |

**The `"semitrans"` driver**: the amplitudes are folded into the occupied-side coefficients first (only one side is transformed into the MO basis — a "semi-transformation"); the transition density factorizes exactly by construction (rank $\le n_\mathrm{occ}$, no SVD needed):

$$
\begin{aligned}
C\!X_{\mu i}^{\mathbb{A}} &= \sum_{a} C_{\mu a}\, z_{ia}^{\mathbb{A}}
&& \text{(eq.1 one batched DGEMM)} \\
K_{\mu\nu}^{\mathbb{A}} &= \sum_{\underline{P},\,i} \bigl(M_{\underline{P}}\, C\!X^{\mathbb{A}}\bigr)_{\mu i}\,\bigl(M_{\underline{P}}\, C_{occ}\bigr)_{\nu i}
&& \text{(eq.2 get\_vk\_ri\_incore\_coeff\_pair)}
\end{aligned}
$$

In eq.2 the right half-transform $M_{\underline{P}} C_{occ}$ is computed once per batch and reused across all trial vectors. Since every $M_{\underline{P}}$ is symmetric, $K[P^{\mathrm{T}}] = K[P]^{\mathrm{T}}$: the B block takes the transpose directly, keeping both sides at $k = n_\mathrm{occ}$. The exchange cost drops from $O(n_\mathrm{aux} n_\mathrm{basis}^3)$ to $O(n_\mathrm{aux} n_\mathrm{basis}^2 n_\mathrm{occ})$. Measured (C6H6/PBE0 TDA, 4 threads): DZ 70.4 → 62.9 s (−11%); TZ 400.1 → 286.6 s (−28%) with ≈230 MB less peak memory; energies agree with `"dm"` to $\sim 10^{-14}$ Ha.

**The `"lowrank"` driver** `get_vk_ri_incore_dm_lowrank` (threshold `tddft_svd_tol`, default 1e-6):

$$
z^{\mathbb{A}} = U \Sigma V^{\mathrm{T}} \qquad
k = \#\{\, \sigma_i \;\ge\; \varepsilon_{\mathrm{svd}}\, \sigma_{\max} \,\}
$$

reducing the exchange cost to $O(n_\mathrm{aux} n_\mathrm{basis}^2 k)$. Note that this is a **lossy approximation**: for typical systems (e.g. benzene/def2-SVP) the trial-vector matrices are effectively full rank, so the low-rank decomposition gives no speedup; it only pays off when transition densities are genuinely low rank. For general use, the exact `"semitrans"` driver is preferred.

## Function `fxc_matvec_ao_batched`

Function path: `ri_tddft::matvec_ao::fxc_matvec_ao_batched`

Applies the XC kernel to the whole trial-vector block in one pass (PySCF `vind(zs)` style). The transition densities are first symmetrized, $P^{\mathrm{sym}} = (P + P^{\mathrm{T}})/2$ (`make_rho_from_dm` assumes a symmetric density for its SIGMA response; the kernel response is invariant under $\mu \leftrightarrow \nu$ symmetrization, so this is exact), then:

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

eq.1 outputs $[n_\mathrm{grid}, n_\mathrm{var}, n_\mathrm{set}]$ and eq.2 outputs $[n_\mathrm{basis}, n_\mathrm{basis}, n_\mathrm{set}]$ (symmetrized internally); both run over all $n_\mathrm{set}$ vectors at once. Coulomb and exchange are flop-bound contractions and remain per-vector calls; the fxc contraction is memory/bandwidth-bound on the grid, and batching means the grid AO values and the kernel table are read only once.

**fxc driver selection** `tddft_fxc_driver` (AO mode only, default `"semitrans"`; `"dm"` is the eq.1–3 path above):

| Value | Dominant per-matvec cost | Extra memory | Remarks |
|--|--|--|--|
| `"semitrans"` (default) | $O(n_\mathrm{occ} n_\mathrm{basis} n_\mathrm{grid})$ (amplitudes folded with $C_{vir}$, contracted against the raw AO) | ψ_occ table only, $(1{+}3\delta_\mathrm{GGA}) n_\mathrm{occ} n_\mathrm{grid}$ | see below |
| `"dm"` | $O(n_\mathrm{basis}^2 n_\mathrm{grid})$ (assembled $[n_\mathrm{basis},n_\mathrm{basis},m]$ transition densities) | none | eq.1–3 path |
| `"mo"` | $O(n_\mathrm{occ} n_\mathrm{vir} n_\mathrm{grid})$ + ψ-table traffic | ψ tables $(n_\mathrm{occ}{+}n_\mathrm{vir})(1{+}3\delta_\mathrm{GGA}) n_\mathrm{grid}$ | see below |

**The `"semitrans"` driver** (default, the `st` branch of `fxc_mo_matvec`): $C_{vir}$ is **folded into the amplitudes** up front so that the virtual side contracts directly against the raw AO values on the grid — no ψ_vir table is ever formed:

$$
\begin{aligned}
\tilde z^{\mathbb{A}}_{i\mu} &= \sum_{a} z^{\mathbb{A}}_{ia}\, C_{\mu a}
&& \text{(eq.S1 one batched DGEMM per call)} \\
\rho_z(g) &= \sum_{i} \psi_i(g) \sum_{\mu} \varphi_\mu(g)\, \tilde z^{\mathbb{A}}_{i\mu}
&& \text{(eq.S2 } \psi_{occ}\text{-vecdot} \times \text{raw-AO GEMM)} \\
v_{1,\alpha}(g) &= w(g)\sum_\beta f^{\mathrm{xc}}_{\alpha\beta}(g)\,\rho_\beta(g)
&& \text{(eq.S3 shared with the other drivers)} \\
E^{\mathbb{A}}_{ia} &= \sum_{\mu} C_{\mu a} \sum_g \varphi_\mu(g)\,\psi_i(g)\, v_{1,\alpha}(g)
&& \text{(eq.S4 two GEMMs: } [\,m n_\mathrm{occ}, n_\mathrm{basis}\,] \times C_{vir}\text{)}
\end{aligned}
$$

It shares the same implementation framework as the `"mo"` driver (`fxc_mo_matvec`: occ-side ψ tables only, grid batches + rayon chunks); the difference is only in the GEMM operands: `"semitrans"` uses the raw AO $[n_\mathrm{batch}, n_\mathrm{basis}]$ as the left operand and the folded amplitudes $[m\, n_\mathrm{occ}, n_\mathrm{basis}]$ as the right one (eq.S1 once per call), plus one extra $[m\, n_\mathrm{occ}, n_\mathrm{basis}] \times C_{vir}$ GEMM on the way back; `"mo"` uses the per-batch projected ψ_vir as the left operand and the raw amplitudes on the right. For UHF the kernel contraction unrolls along the spin-resolved table $f[g,\alpha,\sigma_1,\beta,\sigma_2]$ (the $\sigma$ double loop of eq.S3). Unknown values warn and fall back to `"dm"`.

**The `"mo"` driver**: the AO port of the MO-mode fxc algorithm (mathematically equivalent, implementation-level differences only) — `prepare_ao_data_with_spin` projects and caches the occ/vir MO-on-grid tables in grid batches, and each matrix-vector product contracts purely in the occ/vir space, exactly like MO-mode `fxc_matvec`:

$$\psi_i(g) = \sum_\mu C_{\mu i}\,\varphi_\mu(g) \qquad \psi_a(g) = \sum_\mu C_{\mu a}\,\varphi_\mu(g)$$

(GGA also caches $\partial_d\psi$; layouts are $[n_\mathrm{grid}, \cdot]$). Each matrix-vector product then only contracts in the occ/vir space:

$$
\rho_0(g) = \sum_{ia} z_{ia}\,\psi_i(g)\psi_a(g), \qquad
\rho_{d+1}(g) = \sum_{ia} z_{ia}\,(\partial_d\psi_i\,\psi_a + \psi_i\,\partial_d\psi_a)(g)
$$

$$
v_{1,\alpha}(g) = w(g)\sum_\beta f^{\mathrm{xc}}_{\alpha\beta}(g)\,\rho_\beta(g), \qquad
E_{ia} = \sum_g \Lambda^\alpha_{ia}(g)\,v_{1,\alpha}(g)
$$

with $\Lambda^0_{ia} = \psi_i\psi_a$ and $\Lambda^{d+1}_{ia} = \partial_d\psi_i\,\psi_a + \psi_i\,\partial_d\psi_a$. The output is MO amplitudes directly — no $[n_\mathrm{basis},n_\mathrm{basis},m]$ intermediates, no `contract_back`. (In the current implementation only the occ side is cached; the vir side is streamed per batch — see below.)

Measured trade-off (C6H6/PBE0 TDA, 4 threads):

| | `"dm"` | `"mo"` (chunked) |
|---|---|---|
| DZ | 70.4 s / 909 MB | **46.6 s / 910 MB** |
| TZ | 400.1 s / 2383 MB | **110.7 s / 2150 MB** |

`fxc_mo_matvec` caches only the small occ-side ψ tables (`[ng, nocc]` + gradients); the
vir-side projections are **streamed per grid batch** (one AO evaluation + `C_vir` projection
per batch, batch-local buffers freed afterwards), then the batch runs in rayon-parallel grid
chunks (chunk ≈ 1536, cache-resident `[chunk,·]` buffers; the ρ build contracts z with the
vir side first so the ρ GEMMs run at K = nvir). Energies agree with `"dm"` to ~2e-14 Ha.
`"mo"` is faster at both sizes with comparable-or-lower memory; `"dm"` remains the simpler
fallback (no NIMatmul ψ preparation).

### `grid_batch`: batching over grid points

The full-grid AO tensor $[n_\mathrm{grid}, n_\mathrm{basis}, n_\mathrm{comp}]$ reaches hundreds of MB (or more) even for mid-sized systems. With `grid_batch = true` (default), eq.1–2 run in grid batches (`NIMatmul::split_batch`, the same pattern as the RKS Hessian's `make_hessian_setup_batched`):

$$
\rho^{\mathbb{A}}_{\alpha}(\underline{g}) = \sum_{\mu\nu \in \underline{g}} P^{\mathrm{sym},\mathbb{A}}_{\mu\nu}\, \varphi_\mu(\underline{g})\, \varphi_\nu(\underline{g}) \qquad
F^{\mathbb{A}} \mathrel{{+}{=}} \text{eq.2 contribution on } \underline{g}
$$

The batch size is `NIMatmul.nbatch` (default $1536 \times n_\mathrm{thread}$); the cost is re-evaluating the AO values for each batch at every matrix-vector product (libcint calls are cheap; measured ≈ +15% time for ≈ −40% peak memory; see the developer notes in `tddft_ao.md`).

| Memory type | Expression | Index order | Memory (`grid_batch=false` / `true`) | Remarks |
|--|--|--|--|--|
| fixed | AO value cache | $(g, \mu, c)$ | $n_\mathrm{grid} n_\mathrm{basis} n_\mathrm{comp}$ / $n_\mathrm{batch} n_\mathrm{basis} n_\mathrm{comp}$ | $c$ = density components (1/4) |
| batched | eq.1 output | $(\underline g, \alpha, \mathbb{A})$ | $n_\mathrm{batch} n_\mathrm{var} n_\mathrm{set}$ | |
| batched | eq.2 output | $(\mu, \nu, \mathbb{A})$ | $n_\mathrm{basis}^2 n_\mathrm{set}$ | |
| fixed | kernel table $f^{\mathrm{xc}}_{\alpha\beta}(g)$ | $(g, \alpha, \beta)$ | $n_\mathrm{grid} n_\mathrm{var}^2$ | same in both modes |
| fixed | `out` | $(\mu, \nu, \mathbb{A})$ | $n_\mathrm{basis}^2 n_\mathrm{set}$ | |

## Dense small-system path: `build_a_ao` / `build_b_ao`

Function path: `ri_tddft::matvec_ao::build_a_ao` / `build_b_ao`

When $n_\mathrm{occ} n_\mathrm{vir} \le 15$, `tddft_main` builds the full A (or B) matrix and diagonalizes densely. In AO mode:

$$
A_{ia,jb} = (\varepsilon_a - \varepsilon_i)\delta_{ij}\delta_{ab} + \left[\text{eq.1--3 applied to the identity block } E_{(ia),(jb)}\right]
$$

i.e. the identity matrix is used as the trial-vector block, and one pass through `ao_kernel_block` (folded J + K (including the RSH short-range $K_{SR}$ against `scf.rimatr_sr` with the same K drivers) + batched fxc + contraction back to the MO basis) yields all kernel columns at once, plus the diagonal. The dense path for unrestricted references runs over the concatenated $[\alpha;\beta]$ amplitude space. The mode dispatch is encapsulated in `tddft::build_a`/`build_b` (the MO branch loops unit columns through `a_matvec`), so the solver code is mode-agnostic.

## Batched Davidson interface

AO mode iterates with the batched solvers `solvers::davidson::davidson_solver_batched` / `lr_davidson_solver_batched` (consumer-side alias `tda_davidson_solver_batched`), with closure type:

```rust
FnMut(&MatrixFull<f64>) -> MatrixFull<f64>   // [dim, n_set] → [dim, n_set]
```

The solver hands the whole trial-vector block to the closure, which then: builds the transition densities (one DGEMM covering all $n_\mathrm{set}$ columns) → batched J/K (one `ri_jk` call) → batched fxc (one `make_rho_from_dm` + `make_fxc_pot_with_eff`) → contraction back to the MO basis. Parallelism lives inside the closures (rayon); the subspace iteration itself stays sequential, avoiding nested thread-pool contention.

The per-vector interface `davidson_solver` / `lr_davidson_solver` (alias `tda_davidson_solver`; used by `ri_bse`/`scf_io`) is a per-column adapter around the batched cores.
