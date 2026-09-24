# TDDFT analytic excited-state gradients

Module path: `ri_tddft::tddft_grad` (`src/ri_tddft/tddft_grad.rs`, entry point `TddftGradEngine`)

## Overview

`TddftGradEngine` computes the analytic nuclear gradient of a single excited-state excitation energy. The implementation is a faithful port of PySCF `pyscf/grad/tdrks.py` (`_contract_xc_kernel`, `grad_elec`) and `pyscf/grad/tdrhf.py` (`grad_elec`, the CPHF Z-vector, the Pulay term) to REST's RI representation. The gradient is assembled as

$$
\mathbf{F}_{\mathrm{TDDFT}} = \mathbf{F}_{\mathrm{GS}} + \bigl[\mathbf{F}_{\mathrm{resp}}(x, y) - \mathbf{F}_{\mathrm{resp}}(\mathbf{0}, \mathbf{0})\bigr]
$$

where $\mathbf{F}_{\mathrm{GS}}$ is REST's validated ground-state RKS gradient and the bracketed term is the amplitude-dependent response part (`assemble(x,y) - assemble(0,0)`; the second term cancels the ground-state component of the transition density).

## Theoretical background

With $\omega$ a Rayleigh quotient in the amplitudes (for TDA $y = 0$), the gradient receives three kinds of contributions:

1. **Electronic response gradient**: the transition densities ($\mathbf{P} = C (x+y) C^{\mathrm{T}}$ and difference densities) enter the 0th-order J/K potentials, the derivative J/K integral tables, and two grid passes of the XC kernel (the second pass needs the Z-vector);
2. **CPHF Z-vector**: the relaxed density response is obtained by solving the Z-vector equation (`solve_zvector`, reusing the CPHF machinery of `ri_cphf::CPHFSolverPySCF`);
3. **Pulay terms and nuclear gradients**: additional terms from the position dependence of the auxiliary basis and the grid (`aux_forces_atom`, generator derivatives, overlap derivatives).

The two distinct `get_jk` calls of PySCF are both reproduced in the RI approximation:

- 0th-order `get_jk(mol, dm, hermi=0)`: ordinary J/K AO matrices (`j_potential`/`k_potential`, for the `veff0` series);
- derivative `get_jk`: integral derivative tables $v^{J}_{t,\mu\nu} = -\sum_{ls} \frac{\partial}{\partial R_\mu} (\mu\nu|ls) D_{ls}$ (`vj_bra_atom`/`vk_bra_atom` + `AtomDerivBlocks`, for `veff1`); the auxiliary-basis and metric responses are added per bilinear form from the full potential derivative (`RawRiTensors::d_i_from_blocks` + `d_j_atom`, reusing the infrastructure of `ri_gw::gw_grad`).

## Usage

```toml
[tddft]
     tddft_method =              "tda"    # or "lr"
     nroots =                    3
     tddft_grad_state =          1        # 1-based; 0 disables
```

- Restricted references **only** (panics when `spin_polarization = true`);
- `tddft_grad_state = N` requests the gradient of the N-th excited state; it must not exceed `nroots`, and the TDDFT solve must have run in the same task;
- `main_driver::eval_force` adds the response of the selected state on top of the ground-state gradient and prints the total; geometry optimizations/MD call the same gradient interface, so excited-state structure optimization works directly;
- Channel selection: `tddft_spin = "triplet"` selects the triplet channel; everything else (including `"both"`) uses the singlet channel;
- RSH functionals: the exchange contributions carry the short-range correction (the engine's `hyb` = $c_{LR}$, `hyb_sr` = $c_{SR} - c_{LR}$).

## Key data structures

| Member/struct | Meaning |
|--|--|
| `TddftGradEngine` | per-state gradient engine: `scf`, `state`, `singlet`, `tda`, amplitudes `(x, y)`, raw RI tensors `raw`, pre-inverted RI metric `jinv`, exchange coefficients (`hyb`, `hyb_sr`), fxc cache `fxc_cache`, libcint handle `cint` |
| `RawRiTensors` / `AtomDerivBlocks` | reused from `ri_gw::gw_grad`: atom-derivative blocks of the three-center integrals |
| `FxcHessianCache` / `VindWorkspace` | reused from `dft::response`: the fxc second-order kernel grid cache and evaluation workspace |
| `AOMat` / `Density` / `RankFactor` | column-major AO matrices with an explicit rank factorisation $D = \sum l\, r^{\mathrm{T}}$ (the exchange contractions need the factor structure) |
| `AoBlocks` | optional whole-grid AO table cache (see the memory tuning section) |

## Workflow

```text
TddftGradEngine::new(scf, state, singlet, tda, x, y)
    ├── build raw RI tensors (build_raw_ri_tensors) + metric inverse jinv (once)
    └── response_gradient()
        └── assemble(x, y) - assemble(0, 0)
            ├── transition densities: P = x+y / x-y combinations → dvv, doo, P_pl
            ├── 0th-order potentials: vj / vk (K only when kf ≠ 0) → veff0 series
            ├── XC kernel pass 1: contract_xc_kernel (fxc0)
            ├── Z-vector: solve_zvector(wvo) (CPHF, ri_cphf)
            ├── XC kernel pass 2: contract_xc_kernel (fxcz1, needs the Z-vector)
            ├── derivative potentials: vj_bra_atom / vk_bra_atom (AtomDerivBlocks)
            ├── generator-derivative hcore, overlap derivatives, aux_forces_atom (Pulay)
            └── assemble [3, natm]
```

`contract_xc_kernel` runs twice per gradient (the second pass, `fxcz1`, needs the Z-vector), and both passes walk the same grid with the same `ao_deriv`. `eval_ao_batch` is ~30% of the response cost, so the second pass is pure duplicated work; whether to cache depends on the memory budget (see below).

## Performance and memory tuning (environment variables)

| Environment variable | Purpose |
|--|--|
| `REST_TDDFT_GRAD_AOCACHE_MB` | override the memory budget of the whole-grid AO table cache (cached only when the AO tensor fits; otherwise re-evaluated per batch) |
| `REST_TDDFT_GRAD_NO_AOCACHE` | disable the AO table cache |
| `REST_TDDFT_GRAD_XCBLK_MB` | override the XC grid block size (defaults to the `max_memory` budget) |
| `REST_TDDFT_GRAD_AO_LEGACY=1` | restore the pre-refactor `eval_ao_batch` evaluation path (regression cross-check only; the libcint path is 5-6x faster) |
| `REST_TDDFT_GRAD_TIME` / `REST_TDDFT_GRAD_MEM` | per-step timing / RSS peak-memory tracing on stderr |

## Limitations

- Restricted (RHF/RKS) references only; unrestricted gradients are not implemented;
- Excitation-energy gradients only; no special handling for state crossings (conical intersections);
- The gradient is analytic within the RI approximation (the auxiliary-basis response is included).
