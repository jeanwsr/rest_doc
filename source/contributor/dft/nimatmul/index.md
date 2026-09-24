# nimatmul: DFT Grid Integration and Derivatives

nimatmul (numerical integration by matrix multiplication) is a DFT grid-integration module built on plain matrix multiplications, without einsum or other tensor-contraction tools. The REST implementation lives in `rest/src/dft/numint_matmul` and `rest/src/dft/gen_grids/becke_partitioning_deriv`.

This document set describes the API design concepts, program interface, and design decisions of nimatmul, as well as the conventions and implementation strategy by which it evaluates DFT Hessians (second-order skeleton derivatives of the energy and first-order skeleton derivatives of the Fock matrix), including the grid-shift derivative contributions. The former group of documents was ported from the rstsr-showcase-dft-grids prototype repository, and the Hessian-related documents from the rstsr-showcase-hessian prototype repository; all have been corrected against the official REST implementation.

:::{admonition} AI-assisted translation
:class: note

The English pages in this section were translated from the Chinese originals by an AI agent (Claude Code + glm-5.3 / glm-5.3-flash); all formulas and code are kept verbatim from the Chinese version. The Chinese documents remain the authoritative version.
:::

:::{note}
The `concept`, `basic-api`, and `adr` pages were rewritten and split by an AI agent (Claude Code + glm-5.3) from the design document of rstsr-showcase-dft-grids (cn-dft-api-concepts); the formulas and interfaces have been fact-checked against the current REST implementation, and the notation of the original has been unified to the conventions of `def`.
:::

:::{note}
Unless stated otherwise, each document treats only the closed-shell case; open-shell cases require the corresponding expansion over spin components (see the introductions of each document). Tensor shapes follow the column-major convention; the `becke-deriv` document follows the row-major convention in line with its implementation style, as explained at its beginning.
:::

```{toctree}
:maxdepth: 1
def
concept
basic-api
adr
skeleton2
vmat1
becke-deriv
becke-grid-shift
```
