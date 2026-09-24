# TDDFT Analytic Excited-State Gradients

`TddftGradEngine` (`src/ri_tddft/tddft_grad.rs`) computes the analytic nuclear gradient of a single excited-state excitation energy, triggered by the input keyword `tddft_grad_state`; for user-facing usage see [TDDFT Excited-State Calculations](../../user/tddft.md).

## Implementation notes

- **PySCF reference**: the gradient assembly is a port of `pyscf/grad/tdrks.py` (`_contract_xc_kernel`, `grad_elec`) and `pyscf/grad/tdrhf.py` (`grad_elec`, CPHF Z-vector, Pulay term), rewritten in REST's RI representation.
- **Overall decomposition**: `de_TDDFT = de_GS + response`. `de_GS` reuses REST's validated ground-state RKS gradient; `response = assemble(x, y) - assemble(0, 0)` is the amplitude-dependent part.
- **Two kinds of J/K calls** (as in PySCF): the 0th-order J/K are ordinary AO matrices (corresponding to `veff0doo`/`veff0mop`/`veff0mom`); the derivative tables `vj[t,μ,ν] = -Σ_ls ∂(μν|ls)/∂R D_ls` and `vk[t,μ,ν] = -Σ_ls ∂(μl|νs)/∂R D_ls` (corresponding to `veff1`). Both are reproduced here in the RI approximation.
- **Auxiliary-basis / metric responses**: added per bilinear form from the full potential derivative (`RawRiTensors::d_i_from_blocks` + `d_j_atom`).

> TODO: this page still needs a fuller derivation and numerical-validation notes.
