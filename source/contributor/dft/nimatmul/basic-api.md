# nimatmul Program Interface

This document describes the program interface of the nimatmul module: the interface listing, shape conventions, and key formulas of the pure-function layer and the DFT public-interface layer of the three-layer architecture (see Section 2.4 of [concept.md](concept.md)), as well as the organization of the driver layer.

This document is written against the current REST implementation; the code lives in `rest/src/dft/numint_matmul` and `rest/src/dft/xceff`.

:::{note}
The formula notation and tensor shapes of this document follow the conventions of [def.md](def.md). The list dimension $\mathbb{A}$ is denoted `nset`, the same concept as the arbitrary-property dimension `nprop` of Section 4 in def.md: the length of the density-matrix list, usually corresponding to the number of properties computed, and often standing in for spin $\sigma$ in open-shell density generation.

Shapes applicable to both restricted (spin-unpolarized, closed-shell) and unrestricted (spin-polarized, open-shell) cases are tagged with (R) and (U) respectively for the two cases.

The density-type parameter is named `den_type` in the pure-function and public-interface layers; the same parameter is customarily named `xc_type` in the Hessian driver layer (hess_rks.rs, hess_uks.rs), with identical possible values.
:::

## 1. Module Layout

| File | Architecture layer | Content |
|--|--|--|
| `numint_matmul/pure_eval_rho.rs` | Pure function | Density-on-grids generation `get_rho_from_*_with_output` |
| `numint_matmul/pure_xcpot.rs` | Pure function | XC potential-matrix assembly `{rks,uks}_{vxc,fxc,kxc}_pot_with_eff_with_output` |
| `numint_matmul/nimatmul.rs` | Public interface | `NIMatmul` struct: AO caching, grid batching, density-generation and potential-assembly methods |
| `numint_matmul/hess_rks.rs`, `hess_uks.rs` | Driver | Hessian setup and CP-KS response; see [skeleton2.md](skeleton2.md), [vmat1.md](vmat1.md) |
| `dft/xceff/flags.rs` | Public interface | `XCDenType`, `XCSpin`, `XCPar`, `AO_DERIV_DIM` |
| `dft/xceff/libxc_wrap.rs`, `xc_deriv.rs` | Public interface | `libxc_eval_eff` and the $\gamma \to \rho_r$ expansion of the raw LibXC output |

:::{note}
The `xceff` module is functionally similar to the legacy interface `eval_xc_eff` of `dft/libxc_itrf.rs`; the latter serves `dft/num_int` and `dft/response` (see the interface stack in [../libxc.md](../libxc.md)). The two may be refactored and merged in the future.
:::

## 2. Pure-Function Layer

### 2.1 Density-on-grids generation (`pure_eval_rho.rs`)

The four functions implement the component formulas of Step 1 in Section 2.2 of [concept.md](concept.md); they differ only in the input form: density-matrix inputs perform the GEMM over the basis-function index, bra-ket inputs perform the GEMM over the occupied-orbital index (low-rank optimization). The output is uniformly $\rho_g^{\chi \mathbb{A}}$, i.e. `[ngrids, nvar, nset]`.

Common parameters are `ao` (an AO tensor view of shape `[ngrids, nao, ncomp]`, where `ncomp` is 1/4/4 depending on `den_type`, and 10 for LAPL), `den_type` (`XCDenType`), `out` (pre-allocated output buffer), and `nchunk` (parallel chunk size, see Section 3.1).

**Function `get_rho_from_dm_with_output`**

Density-matrix input. For RHO/SIGMA only one GEMM is needed (eq.1); for TAU one additional GEMM per spatial component $r$ (the $\bar{\phi}_{g \mu}^{(r)}$ of eq.4). **The density matrices must be symmetric**: the GEMM performs only a one-sided contraction, with the other side done by the scaled reduction; symmetry guarantees that this is equivalent to the full quadratic form.

$$
\begin{aligned}
\bar{\phi}_{g \mu} &= \sum_\nu D_{\mu \nu}^{\mathbb{A}} \, \phi_{g \nu} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_\mu \phi_{g \mu} \, \bar{\phi}_{g \mu} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= 2 \sum_\mu \phi_{g \mu}^r \, \bar{\phi}_{g \mu} && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_\mu \phi_{g \mu}^r \, \bar{\phi}_{g \mu}^{(r)} && \text{(eq.4)}
\end{aligned}
$$

where $\bar{\phi}_{g \mu}^{(r)} = \sum_\nu D_{\mu \nu}^{\mathbb{A}} \, \phi_{g \nu}^r$.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `dm_list` | $D_{\mu \nu}^{\mathbb{A}}$ | $(\mu, \nu)$<br>`[u, v]` | `nset` matrices of `[nao, nao]` |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

| Memory type | Equation | Expression | Index order | Size |
|--|--|--|--|--|
| thread | (eq.1) | $\bar{\phi}_{g \mu}$ | $(g, \mu)$ | `(nchunk, nao)` |
| thread | (eq.2)–(eq.4) | local output | $(g, \chi)$ | `(nchunk, nvar)` |

The $\bar{\phi}_{g \mu}$ of eq.1 is reused by eq.2/eq.3; for TAU it is recomputed for each $r$. Peak memory is `nthreads × nchunk × (nao + nvar)`.

**Function `get_rho_from_homogeneous_braket_with_output`**

Homogeneous bra-ket input, i.e. $D^{\mathbb{A}} = C^{\mathbb{A}} (C^{\mathbb{A}})^T$. The occupation information should be folded into the coefficients: for a closed shell with occupation 2, the user should pass $C_{\mu i} \sqrt{2}$ (or multiply the output by 2 afterwards; we prefer the former, and the REST driver indeed folds $\sqrt{n_i}$ into the occupied-orbital coefficients).

$$
\begin{aligned}
\phi_{g i}^{\mathbb{A}} &= \sum_\mu \phi_{g \mu} \, C_{\mu i}^{\mathbb{A}} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_i \phi_{g i}^{\mathbb{A}} \, \phi_{g i}^{\mathbb{A}} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= 2 \sum_i \phi_{g i}^{r, \mathbb{A}} \, \phi_{g i}^{\mathbb{A}} && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_i \phi_{g i}^{r, \mathbb{A}} \, \phi_{g i}^{r, \mathbb{A}} && \text{(eq.4)}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra_list` | $C_{\mu i}^{\mathbb{A}}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- The `nocc` of each set in the list may differ.

| Memory type | Equation | Expression | Index order | Size |
|--|--|--|--|--|
| thread | (eq.1) | orbital-on-grids buffers ×2 | $(g, i)$ | `(nchunk, nocc_max)` |
| thread | (eq.2)–(eq.4) | local output | $(g, \chi)$ | `(nchunk, nvar)` |

Peak memory is `nthreads × nchunk × (2 nocc_max + nvar)`.

**Function `get_rho_from_one_bra_mult_ket_with_output`**

Bra-ket input with one shared bra and multiple kets. This is the typical form of the $X_{ai}^{\mathbb{A}}$-type perturbations in TD/CP-KS equations: for the perturbation density $D_{\mu \nu}^{\mathbb{A}} = \sum_{ai} C_{\mu i}^{\text{bra}} X_{ai}^{\mathbb{A}} C_{\nu a}^{\text{ket}}$, the bra side can be passed as the occupied-orbital coefficients, while the ket side is best half-transformed, $\tilde{C}_{\nu i}^{\mathbb{A}, \text{ket}} = \sum_a X_{ai}^{\mathbb{A}} C_{\nu a}^{\text{ket}}$, before being passed in.

$$
\begin{aligned}
\phi_{g i}^{(\text{bra})} &= \sum_\mu \phi_{g \mu} \, C_{\mu i}^{(\text{bra})}, \quad \phi_{g i}^{(\text{ket}, \mathbb{A})} = \sum_\mu \phi_{g \mu} \, C_{\mu i}^{(\text{ket}, \mathbb{A})} && \text{(eq.1)} \\
\xi_g^{\chi = \rho, \, \mathbb{A}} &= \sum_i \phi_{g i}^{(\text{bra})} \phi_{g i}^{(\text{ket}, \mathbb{A})} && \text{(eq.2)} \\
\xi_g^{\chi = \rho_r, \, \mathbb{A}} &= \sum_i \left( \phi_{g i}^{r, (\text{bra})} \phi_{g i}^{(\text{ket}, \mathbb{A})} + \phi_{g i}^{(\text{bra})} \phi_{g i}^{r, (\text{ket}, \mathbb{A})} \right) && \text{(eq.3)} \\
\xi_g^{\chi = \tau, \, \mathbb{A}} &= \sum_r \frac{1}{2} \sum_i \phi_{g i}^{r, (\text{bra})} \phi_{g i}^{r, (\text{ket}, \mathbb{A})} && \text{(eq.4)}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra` | $C_{\mu i}^{(\text{bra})}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{(\text{ket}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- The `nocc` of every matrix in `ket_list` must equal that of `bra`.

Peak memory is `nthreads × nchunk × (3 nocc + nvar)`.

**Function `get_rho_from_mult_bra_mult_ket_with_output`**

Bra-ket input with multiple bras and multiple kets. The formulas are identical to `get_rho_from_one_bra_mult_ket_with_output`; only the bra now varies with the set $\mathbb{A}$ as well.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `ao` | $\phi_{g \mu}^{*}$ | $(g, \mu, *)$<br>`[g, u, *]` | `[ngrids, nao, ncomp]` |
| `bra_list` | $C_{\mu i}^{(\text{bra}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{(\text{ket}, \mathbb{A})}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| `out`<br>(output) | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `nchunk` | | | |

- `bra_list` and `ket_list` have the same length, with one-to-one corresponding occupations.

Peak memory is `nthreads × nchunk × (3 nocc_max + nvar)`.

### 2.2 XC potential-matrix assembly (`pure_xcpot.rs`)

All functions share two private contraction cores, i.e. the `contract_ao_wv` of Step 3 in Section 2.2 of [concept.md](concept.md):

- `contract_ao_wv_without_symmetrize` produces the asymmetric half matrix (terms selected by `den_type`):

    $$
    \begin{aligned}
    (\text{half})_{\mu \nu} = \;& 0.5 \sum_g \phi_{g \mu} \, w_g f_g^{\rho} \, \phi_{g \nu} \\
    &+ \sum_{t} \sum_g \phi_{g \mu}^t \, w_g f_g^{\rho_t} \, \phi_{g \nu} \\
    &+ 0.25 \sum_{t} \sum_g \phi_{g \mu}^t \, w_g f_g^{\tau} \, \phi_{g \nu}^t
    \end{aligned}
    $$

    followed by the uniform symmetrization $V \leftarrow (\text{half}) + (\text{half})^T$. The reasons for the coefficients 0.5/1.0/0.25 are given in Section 3 of [adr.md](adr.md).

- `contract_ao_wv_bra` takes the bra-transformed AO and directly produces the asymmetric `[nao, nocc]` output without symmetrization, with coefficients 1.0/1.0/0.5 respectively (see eq.3 of the bra-trans function below).

The public functions differ only in how the "effective potential" $wv$ entering the contraction core is composed. Common parameters are `den_type`, `ao`, `weights` ($w_g$, `[ngrids]`), and `nchunk`.

**Function `rks_vxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi} &= w_g f_g^{\chi} && \text{(eq.1)} \\
V_{\mu\nu}^{\text{xc}} &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{(eq.2)}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `vxc_eff` | $f_g^{\chi}$ | $(g, \chi)$<br>`[g, x]` | `[ngrids, nvar]` |
| `vxc`<br>(output) | $V_{\mu\nu}^{\text{xc}}$ | $(\mu, \nu)$<br>`[u, v]` | `[nao, nao]` |

- `vxc_eff` is the functional-evaluation output; the density used is generally the self-consistent-field density.
- The output $V_{\mu\nu}^{\text{xc}}$ is symmetric.

| Memory type | Expression | Index order | Size |
|--|--|--|--|
| thread | `wv` contraction buffer | $(g, \mu)$ | `(nchunk, nao)` |
| thread | local output | $(\mu, \nu)$ | `(nao, nao)` |

The per-thread local outputs are accumulated into the global output under a lock.

**Function `rks_fxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi, \mathbb{A}} &= w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_g^{\chi'}[\mathbf{R}^{\mathbb{A}}] && \text{(eq.1)} \\
F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}] &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{(eq.2)}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `fxc`<br>(output) | $F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \mathbb{A})$<br>`[u, v, A]` | `[nao, nao, nset]` |

- The output is symmetric for each $\mathbb{A}$.

Memory usage is the same as `rks_vxc_pot_with_eff_with_output`.

**Function `rks_kxc_pot_with_eff_with_output`**

$$
\begin{aligned}
wv_g^{\chi, (\mathbb{A}, \mathbb{B})} &= w_g \sum_{\chi' \chi''} f_g^{\chi\chi'\chi''} \, \xi_g^{\chi'}[\mathbf{R}'^{\mathbb{A}}] \, \xi_g^{\chi''}[\mathbf{R}''^{\mathbb{B}}] && \text{(eq.1)} \\
K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}] &= (\text{half})_{\mu\nu} + (\text{half})_{\nu\mu} && \text{eq.2}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \chi' \chi''}$ | $(g, \chi, \chi', \chi'')$<br>`[g, x, y, z]` | `[ngrids, nvar, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset1]` |
| `rho2` | $\xi_g^{\chi}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \mathbb{B})$<br>`[g, x, B]` | `[ngrids, nvar, nset2]` |
| `kxc`<br>(output) | $K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \mathbb{A}, \mathbb{B})$<br>`[u, v, A, B]` | `[nao, nao, nset1, nset2]` |

- The list lengths `nset1` and `nset2` of `rho1` and `rho2` may differ.
- The output is symmetric for each $(\mathbb{A}, \mathbb{B})$.

**Function `rks_fxc_pot_with_eff_bra_trans_with_output`**

The bra-transformed fxc: the AO are first contracted with the occupied-orbital coefficients, producing the half-transformed `[nao, nocc]` output; no symmetrization is needed.

$$
\begin{aligned}
\phi_{g i}^{(\text{bra}), *} &= \sum_\mu \phi_{g \mu}^{*} \, C_{\mu i} \quad (\text{for each component } *) && \text{(eq.1)} \\
wv_g^{\chi, \mathbb{A}} &= w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_g^{\chi'}[\mathbf{R}^{\mathbb{A}}] && \text{(eq.2)} \\
F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]
&= \sum_g \phi_{g \mu} \, wv_g^{\rho, \mathbb{A}} \, \phi_{g i}^{(\text{bra})} \\
&\quad + \sum_{t,g} \left( \phi_{g \mu}^t \, wv_g^{\rho_t, \mathbb{A}} \, \phi_{g i}^{(\text{bra})} + \phi_{g \mu} \, wv_g^{\rho_t, \mathbb{A}} \, \phi_{g i}^{t, (\text{bra})} \right) \\
&\quad + \frac{1}{2} \sum_{t,g} \phi_{g \mu}^t \, wv_g^{\tau, \mathbb{A}} \, \phi_{g i}^{t, (\text{bra})} && \text{(eq.3)}
\end{aligned}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `bra` | $C_{\mu i}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `fxc`<br>(output) | $F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$<br>`[u, i, A]` | `[nao, nocc, nset]` |

- The `nset` of `rho1` is the number of kets.

| Memory type | Expression | Index order | Size |
|--|--|--|--|
| persistent | $\phi_{g i}^{(\text{bra}), *}$ (eq.1) | $(g, i, *)$ | `(ngrids, nocc, ncomp)` |
| thread | `wv` contraction buffer | $(g, i)$ | `(nchunk, nocc)` |
| thread | local output | $(\mu, i)$ | `(nao, nocc)` |

The bra orbital-on-grids of eq.1 is generated once inside the function and held throughout.

**Function `uks_vxc_pot_with_eff_with_output`**

The formulas are identical to `rks_vxc_pot_with_eff_with_output` (eq.1/eq.2) with $f_g^{\chi} \to f_g^{\chi \sigma}$; the contraction and symmetrization are performed independently for each spin channel $\sigma$.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `vxc_eff` | $f_g^{\chi \sigma}$ | $(g, \chi, \sigma)$<br>`[g, x, σ]` | `[ngrids, nvar, 2]` |
| `vxc`<br>(output) | $V_{\mu\nu}^{\text{xc}, \sigma}$ | $(\mu, \nu, \sigma)$<br>`[u, v, σ]` | `[nao, nao, 2]` |

**Function `uks_fxc_pot_with_eff_with_output`**

The formulas are identical to `rks_fxc_pot_with_eff_with_output`, with the grid-space contraction including the spin channel $\sum_{\chi' \sigma'} f^{\chi\sigma, \chi'\sigma'} \xi^{\chi'\sigma'}[\mathbf{R}^{\mathbb{A}}]$.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \sigma \chi' \sigma'}$ | $(g, \chi, \sigma, \chi', \sigma')$<br>`[g, x, σ, y, ς]` | `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$<br>`[g, x, σ, A]` | `[ngrids, nvar, 2, nset]` |
| `fxc`<br>(output) | $F_{\mu\nu}^{\text{xc}, \sigma} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \sigma, \mathbb{A})$<br>`[u, v, σ, A]` | `[nao, nao, 2, nset]` |

- The output is symmetric for each $(\sigma, \mathbb{A})$.

**Function `uks_kxc_pot_with_eff_with_output`**

The formulas are identical to `rks_kxc_pot_with_eff_with_output`, with a spin dimension added to each tensor.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \sigma \chi' \sigma' \chi'' \sigma''}$ | $(g, \chi, \sigma, \chi', \sigma', \chi'', \sigma'')$ | `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$ | `[ngrids, nvar, 2, nset1]` |
| `rho2` | $\xi_g^{\chi \sigma}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \sigma, \mathbb{B})$ | `[ngrids, nvar, 2, nset2]` |
| `kxc`<br>(output) | $K_{\mu\nu}^{\text{xc}, \sigma} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \sigma, \mathbb{A}, \mathbb{B})$ | `[nao, nao, 2, nset1, nset2]` |

**Function `uks_fxc_pot_with_eff_bra_trans_with_output`**

The formulas are identical to `rks_fxc_pot_with_eff_bra_trans_with_output`, but the bra and the output are two independent tensors for the $\alpha, \beta$ spins; since the input and output types differ, it cannot be merged with the closed-shell function (see Section 3.4).

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \sigma \chi' \sigma'}$ | $(g, \chi, \sigma, \chi', \sigma')$ | `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi \sigma}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \sigma, \mathbb{A})$ | `[ngrids, nvar, 2, nset]` |
| `bra` | $C_{\mu i}^{(\text{bra}, \sigma)}$ | $(\mu, i)$ | two tensors `[nao, nocc_α]`, `[nao, nocc_β]` |
| `fxc`<br>(output) | $F_{\mu i}^{\text{xc}, \sigma} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$ | two tensors `[nao, nocc_σ, nset]` |

## 3. DFT Public-Interface Layer

### 3.1 The `NIMatmul` Struct

`NIMatmul` is the main grid-driven struct, defined in `nimatmul.rs`. Its key fields:

| Field | Description |
|--|--|
| `cint` | Integral engine (libcint wrapper) |
| `coords`, `weights` | Grid coordinates and weights |
| `atm_idx` | Atom index each grid point belongs to; consumed by the Becke grid-shift derivatives (`usize::MAX` for grids attached to no atom) |
| `quadrature_weights` | (radial × angular) quadrature weights before the Becke partitioning; likewise consumed by the grid-shift derivatives |
| `cache_tensor` | AO cache, keyed by derivative order (e.g. `"ao_deriv0"`, `"ao_deriv1"`), holding copy-on-write tensors `TsrCow` |
| `nchunk` | Parallel chunk size; default 1536, usually a few times the KC dimension of GEMM micro-kernels |
| `nbatch` | Memory batch size; default $1536 \times 1 \times n_\text{threads}$ |

`NIMatmul::new` takes five inputs: `cint`, `coords`, `weights`, `atm_idx`, and `quadrature_weights`; the last two are prepared for the Becke grid-shift derivatives (see [becke-grid-shift.md](becke-grid-shift.md)). **The grid-shift derivatives require atom-grouped grids**; if the grids were generated by multiple threads with interleaved atom indices, regroup them with `regroup_grids_by_atom` first (see [becke-grid-shift.md](becke-grid-shift.md)).

**Grid batching vs grid chunking**: `nbatch` controls memory usage (the full AO tensor `[ngrids, nao, ncomp]` can be too large), while `nchunk` controls parallel granularity. The relation is full-grid > batch > chunk > per-grid = 1.

**AO cache strategy**: `get_cached_ao(deriv)` computes and caches the AO when needed. If a higher derivative order is already cached, the lower-order subset can be sliced out of it without recomputation (this is also why the cache holds `TsrCow`, i.e. copy-on-write tensors). `prepare_ao(deriv)` computes the AO on grids through libcint's `eval_gto`, with output shape `[ngrids, nao, ncomp]`. `split_batch(start, end)` slices a batch instance out of the full-grid instance, slicing the cached AO tensors accordingly.

### 3.2 Density-on-grids Generation Methods

| Method | Description | Typical scenario |
|--|--|--|
| `make_rho_from_dm` | Generate the density from a list of density matrices $D_{\mu \nu}^{\mathbb{A}}$ | Basic functionality; relaxed-density computations in post-SCF |
| `make_rho_from_homogeneous_braket` | Generate the density from homogeneous bra-ket coefficients $C_{\mu i}^{\mathbb{A}}$ | Self-consistent field |
| `make_rho_from_one_bra_mult_ket` | Shared bra $C_{\mu i}^{\text{bra}}$, multiple kets $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | Gradient properties |
| `make_rho_from_mult_bra_mult_ket` | Multiple bras $C_{\mu i}^{\mathbb{A}, \text{bra}}$, multiple kets $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | |

The four methods wrap the same-named pure functions of Section 2.1: `&mut self` internally handles the AO cache (`get_cached_ao`) and the choice of `nchunk`, and the output is a newly allocated tensor. Formulas, input requirements, and shapes are given in Section 2.1; here we only add API usage notes.

**Function `make_rho_from_dm`**

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `dm_list` | $D_{\mu \nu}^{\mathbb{A}}$ | $(\mu, \nu)$<br>`[u, v]` | `nset` matrices of `[nao, nao]` |
| `den_type` | | | `XCDenType` |
| Return value | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- For closed-shell self-consistent-field computations, a single density must be converted into a length-1 list.
- For open-shell computations, the property index $\mathbb{A}$ can stand in for the spin index $\sigma$: passing the list `[dm_α, dm_β]` yields the output `[ngrids, nvar, 2]`, exactly the open-shell density shape.

**Function `make_rho_from_homogeneous_braket`**

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `bra_list` | $C_{\mu i}^{\mathbb{A}}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| Return value | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- The occupations of the sets in the list need not be identical: one may pass coefficient matrices with different $\alpha, \beta$ electron counts for the open shell.
- The closed-shell occupation is usually 2; the user should either multiply $C_{\mu i}$ by $\sqrt{2}$ before passing them in, or multiply the output by 2 afterwards. We prefer the former.

**Function `make_rho_from_one_bra_mult_ket`**

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `bra` | $C_{\mu i}^{\text{bra}}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| Return value | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- This function is motivated by property computations: when solving the TD/CP-KS equations, the $X_{ai}^{\mathbb{A}}$ involved must first be transformed into the atomic-orbital representation before entering the DFT computation; since occupied orbitals are usually far fewer than AOs, passing $C_{\mu i}^{\text{bra}}$ saves work, while the other side should be half-transformed (see the same-named pure function in Section 2.1).
- For the open shell, this function cannot handle $\alpha$ and $\beta$ spins simultaneously; two separate calls are needed.

**Function `make_rho_from_mult_bra_mult_ket`**

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `bra_list` | $C_{\mu i}^{\mathbb{A}, \text{bra}}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `ket_list` | $C_{\mu i}^{\mathbb{A}, \text{ket}}$ | $(\mu, i)$<br>`[u, i]` | `nset` matrices of `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| Return value | $\xi_g^{\chi \mathbb{A}}$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |

- `bra_list` and `ket_list` must have the same length, with one-to-one corresponding occupations.

### 3.3 Functional-Evaluation Bridge (`xceff`)

`dft/xceff/flags.rs` defines three enums and one constant:

| Name | Values | Description |
|--|--|--|
| `XCDenType` | `RHO` / `SIGMA` / `TAU` / `LAPL` | Density type; progressive design, see Section 7 of [adr.md](adr.md) |
| `XCSpin` | `Unpolarized` / `Polarized` | Spin-polarized or not |
| `XCPar` | `Par { chunk_size }` / `Serial` | Parallelization strategy of functional evaluation; convertible from `usize`, `Option<usize>`, `bool` |
| `AO_DERIV_DIM` | `[1, 4, 10, 20, 35]` | AO component counts of each derivative order (up to fourth) |

`determine_den_type` decides the density type from the functional family: LDA/HybLDA → `RHO`, GGA/HybGGA → `SIGMA`, mGGA/HybMGGA → `TAU` or `LAPL` depending on `needs_laplacian()`; `determine_den_type_from_list` takes the strictest type across a functional list.

**Function `libxc_eval_eff`**

Converts the raw LibXC output into the effective-potential format. LibXC uses $\gamma$ as its GGA variable; the chain-rule expansion $\gamma \to \rho_r$ (including the diagonal correction terms of second order and above) happens inside the function via `transform_xc_inner` of `xc_deriv.rs` (see Section 5 of [adr.md](adr.md)).

```rust
pub fn libxc_eval_eff(
    xc_func: &LibXCFunctional,
    rho: TsrView,
    deriv: usize,
    par: impl Into<XCPar>,
) -> Vec<Tsr>
```

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `xc_func` | | | `LibXCFunctional` |
| `rho` | (R) $\xi_g^{\chi}$ <br> (U) $\xi_g^{\chi \sigma}$ | $(g, \chi)$<br>`[g, x]` | (R) `[ngrids, nvar]` <br> (U) `[ngrids, nvar, 2]` |
| `deriv` | | | |
| `par` | | | `XCPar` |
| Return value | $f_g$ and derivatives of all orders | | `Vec<Tsr>` |

- The `LibXCFunctional` instance must be created with its spin specified; the program uses the instance to decide between closed/open shell.
- `deriv` is the **highest** derivative order to compute; the return value contains the outputs of all orders from 0 to `deriv`, the $k$-th element being the order-$k$ output.

The effective-potential output shapes are:

| `deriv` | Closed shell | Open shell |
|--|--|--|
| 0 | `[ngrids]` (exc, $f_g$) | `[ngrids]` |
| 1 | `[ngrids, nvar]` (vxc_eff, $f_g^{\chi}$) | `[ngrids, nvar, 2]` |
| 2 | `[ngrids, nvar, nvar]` (fxc_eff, $f_g^{\chi \chi'}$) | `[ngrids, nvar, 2, nvar, 2]` |
| 3 | `[ngrids, nvar, nvar, nvar]` (kxc_eff, $f_g^{\chi \chi' \chi''}$) | `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |

The default parallel chunk sizes depend on the density type and spin: RHO/unpolarized 16384, RHO/polarized 6144, SIGMA 384, TAU/LAPL 256. **Nested parallelism degrades to serial**: if the calling thread is already inside a rayon thread pool, functional evaluation uses a single thread, to avoid thread oversubscription.

### 3.4 XC Potential-Matrix Assembly Methods

| Method | Description |
|--|--|
| `make_vxc_pot_with_eff` | First-order XC potential |
| `make_fxc_pot_with_eff` | Second-order XC kernel |
| `make_kxc_pot_with_eff` | Third-order XC kernel |
| `make_rks_fxc_pot_with_eff_bra_trans` | Second-order XC kernel (bra-transformed, low-rank optimized) |
| `make_uks_fxc_pot_with_eff_bra_trans` | Second-order XC kernel (bra-transformed, low-rank optimized, open shell) |

The five methods wrap the same-named pure functions of Section 2.2: `&mut self` internally handles the AO cache and the multiplication of the grid weights $w_g$, so $w_g$ and the AO are not passed as parameters. Formulas and input requirements are given in Section 2.2; here we present the unified formulas and API usage notes.

**Function `make_vxc_pot_with_eff`**

$$
V_{\mu \nu}^\text{xc} = \sum_g \sum_\chi w_g f_g^\chi \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `vxc_eff` | (R) $f_g^{\chi}$ <br> (U) $f_g^{\chi \sigma}$ | $(g, \chi)$<br>`[g, x]` | (R) `[ngrids, nvar]` <br> (U) `[ngrids, nvar, 2]` |
| `den_type` | | | `XCDenType` |
| `spin` | | | `XCSpin` |
| Return value | $V_{\mu\nu}^{\text{xc}}$ | $(\mu, \nu)$<br>`[u, v]` | (R) `[nao, nao]` <br> (U) `[nao, nao, 2]` |

- `spin` decides whether the RKS or UKS pure function is called.

$\partial \xi_g^{\chi} / \partial D_{\mu\nu}$ is independent of the density matrix; it depends only on the orbital-on-grids $\phi_{g \mu}$ and their gradients. Like $w_g$, it is managed through `&mut self`.

**Function `make_fxc_pot_with_eff`**

$$
F_{\mu\nu}^\text{xc} [\mathbf{R}^{\mathbb{A}}] = \sum_g \sum_\chi \left( w_g \sum_{\chi'} f_g^{\chi\chi'} \, \xi_{g}^{\chi'}[\mathbf{R}^{\mathbb{A}}] \right) \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | (R) `[ngrids, nvar, nvar]` <br> (U) `[ngrids, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | (R) `[ngrids, nvar, nset]` <br> (U) `[ngrids, nvar, 2, nset]` |
| `den_type`, `spin` | | | |
| Return value | $F_{\mu\nu}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, \nu, \mathbb{A})$<br>`[u, v, A]` | (R) `[nao, nao, nset]` <br> (U) `[nao, nao, 2, nset]` |

**Function `make_kxc_pot_with_eff`**

$$
K_{\mu\nu}^\text{xc} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}] = \sum_g \sum_\chi \left( w_g \sum_{\chi' \chi''} f_g^{\chi\chi'\chi''} \, \xi_{g}^{\chi'}[\mathbf{R}'^{\mathbb{A}}] \, \xi_{g}^{\chi''}[\mathbf{R}''^{\mathbb{B}}] \right) \frac{\partial \xi_g^{\chi}}{\partial D_{\mu\nu}}
$$

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `kxc_eff` | $f_g^{\chi \chi' \chi''}$ | $(g, \chi, \chi', \chi'')$<br>`[g, x, y, z]` | (R) `[ngrids, nvar, nvar, nvar]` <br> (U) `[ngrids, nvar, 2, nvar, 2, nvar, 2]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}'^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | (R) `[ngrids, nvar, nset1]` <br> (U) `[ngrids, nvar, 2, nset1]` |
| `rho2` | $\xi_g^{\chi}[\mathbf{R}''^{\mathbb{B}}]$ | $(g, \chi, \mathbb{B})$<br>`[g, x, B]` | (R) `[ngrids, nvar, nset2]` <br> (U) `[ngrids, nvar, 2, nset2]` |
| `den_type`, `spin` | | | |
| Return value | $K_{\mu\nu}^{\text{xc}} [\mathbf{R}'^{\mathbb{A}}, \mathbf{R}''^{\mathbb{B}}]$ | $(\mu, \nu, \mathbb{A}, \mathbb{B})$<br>`[u, v, A, B]` | (R) `[nao, nao, nset1, nset2]` <br> (U) `[nao, nao, 2, nset1, nset2]` |

**Function `make_rks_fxc_pot_with_eff_bra_trans`**

This function effectively computes

$$
F_{\mu i}^\text{xc} [\mathbf{R}^{\mathbb{A}}] = \sum_\nu F_{\mu\nu}^\text{xc} [\mathbf{R}^{\mathbb{A}}] C_{\nu i}
$$

In practice, however, the best-performing pattern is generally to low-rank-optimize the product of $\partial \xi_g^{\chi} / \partial D_{\mu\nu}$ with $C_{\nu i}$ (the `contract_ao_wv_bra` of Section 2.2).

The reason such a function is needed is that the TD/CP-KS equations frequently involve computations of the form

$$
A_{a i, b j}^\text{xc} R_{b j} = \sum_{\mu \nu} C_{\mu a} C_{\nu i} F_{\mu\nu}^\text{xc} [\mathbf{R}]
$$

and contracting the occupied orbitals first to obtain $F_{\mu i}^\text{xc} [\mathbf{R}]$ performs better, so the above becomes

$$
A_{a i, b j}^\text{xc} R_{b j} = \sum_{\mu i} C_{\mu a} F_{\mu i}^\text{xc} [\mathbf{R}]
$$

Obtaining $F_{\mu i}^\text{xc} [\mathbf{R}]$ is only a half-transformation. Note that the full transformation barely reduces the cost at $n_\mathrm{set} = 1$, and increases the cost for $n_\mathrm{set} > 1$; hence for CP-KS situations (multiple property matrices to compute), the full transformation is inferior to the half one. This is why only the half-transformed function is provided at the public-interface layer.

| Variable | Meaning | Index order | Shape |
|--|--|--|--|
| `fxc_eff` | $f_g^{\chi \chi'}$ | $(g, \chi, \chi')$<br>`[g, x, y]` | `[ngrids, nvar, nvar]` |
| `rho1` | $\xi_g^{\chi}[\mathbf{R}^{\mathbb{A}}]$ | $(g, \chi, \mathbb{A})$<br>`[g, x, A]` | `[ngrids, nvar, nset]` |
| `bra` | $C_{\mu i}$ | $(\mu, i)$<br>`[u, i]` | `[nao, nocc]` |
| `den_type` | | | `XCDenType` |
| Return value | $F_{\mu i}^{\text{xc}} [\mathbf{R}^{\mathbb{A}}]$ | $(\mu, i, \mathbb{A})$<br>`[u, i, A]` | `[nao, nocc, nset]` |

- The `nset` of `rho1` is the number of kets in the bra-ket form.
- **This function is for the closed shell only**. The reason is that the open-shell input `bra` and output cannot be a single higher-order tensor carrying a spin dimension; they must be two separate $\alpha$ and $\beta$ tensors. Since the types differ, a single function cannot realize both in strongly typed languages such as Rust.

**Function `make_uks_fxc_pot_with_eff_bra_trans`**

Same meaning as above, but for the open shell only. The `bra` are the two spin coefficient matrices `[nao, nocc_α]` and `[nao, nocc_β]`, and the return value is the pair `[nao, nocc_σ, nset]` for the two spins. Since the input and output types differ, it cannot be merged with the closed-shell function.

## 4. Driver Layer

In REST, the XC computations of SCF and gradients are currently handled by `dft/num_int` (with `libxc_itrf`); the driver layer of nimatmul targets the **DFT Hessian and its CP-KS response**, living in `hess_rks.rs` and `hess_uks.rs`. The formula derivations and implementation strategies are given in [skeleton2.md](skeleton2.md) (second-order skeleton derivatives of the energy) and [vmat1.md](vmat1.md) (first-order skeleton derivatives of the Fock matrix); here we only list the interface organization.

| Function / struct | Description |
|--|--|
| `get_rho_exc_vxc_fxc` (`_uks`) | From the AO and the occupied-contracted AO-density matrices, obtain `rho`, `exc`, `vxc`, `fxc` |
| `eval_vxc_fxc_from_rho` (`_uks`) | From an existing `rho`, compute `vxc_eff` and `fxc_eff` directly (skipping density generation) |
| `make_cpks_vxc_fxc` (`_uks`) | Lean vxc/fxc evaluation: density from bra-ket, minimal required AO derivative order, for a dedicated CP-KS grid |
| `make_hessian_setup_becke` (`_uks`) | Main Hessian setup: density generation, functional evaluation, and the skeleton intermediates inside a grid-batched loop |
| `get_rks_response_bra` / `get_rks_response_bra_batched` (`_uks`) | CP-KS response: `make_rho_from_one_bra_mult_ket` + `make_fxc_pot_with_eff_bra_trans` |
| `RHessKSNIMatmul` / `UHessKSNIMatmul` | Driver structs; hold the `NIMatmul`, the functional list, and intermediates, plugging into the Hessian solver through traits such as `HessUtilAPI`, `RHessElecInteractAPI` |

The common pattern of the driver layer is to batch the grids by `nbatch` and chain the three steps within each batch (the CP-KS response as an example, `get_rks_response_bra_batched`):

```rust
for start in (0..ngrids).step_by(nbatch) {
    let end = (start + nbatch).min(ngrids);
    // 1. slice out the batched NIMatmul (cached AO tensors are sliced accordingly)
    let mut ni_batch = ni.split_batch(start, end);
    // 2. density-on-grids generation
    let rho1 = ni_batch.make_rho_from_one_bra_mult_ket(mocc, &mo1_bra_list, den_type);
    // 3. XC potential-matrix assembly (fxc_eff comes from functional evaluation, sliced per batch)
    let resp_batch = ni_batch.make_rks_fxc_pot_with_eff_bra_trans(
        fxc_eff.i(start..end), rho1, mocc, den_type);
    resp += resp_batch;
}
```

Functional evaluation (`libxc_eval_eff`) happens once on the full grid, outside the batch loop; the batch loop only slices the AO, generates densities, and assembles potentials. For non-response computations such as vxc, the accumulation of the energy and electron count ($E^\text{xc} = \sum_g w_g f_g \rho_g$, $N_e = \sum_g w_g \rho_g$) also happens in this layer.

:::{note}
In the closed-shell CP-KS response, the occupied-orbital coefficients are not multiplied by $\sqrt{2}$; the occupation contribution is restored by a factor of 4 on the output (`4.0 * resp` in `get_rks_response_bra`). This is the same matter as the $\sqrt{n_i}$ strategy of `make_rho_from_homogeneous_braket` in Section 3.2, implemented differently.
:::
