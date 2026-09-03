# nimatmul: DFT Grid Integration and Derivatives

nimatmul (numerical integration by matrix multiplication) is a DFT grid-integration module built on plain matrix multiplications, without einsum or other tensor-contraction tools. The REST implementation lives in `rest/src/dft/numint_matmul` and `rest/src/dft/gen_grids/becke_partitioning_deriv`.

This document set describes the conventions and implementation strategy by which nimatmul evaluates DFT Hessians (second-order skeleton derivatives of the energy and first-order skeleton derivatives of the Fock matrix), including the grid-shift derivative contributions. The documents were ported from the rstsr-showcase-hessian prototype repository and corrected against the official REST implementation.

:::{admonition} AI-assisted translation
:class: note

The English pages in this section were translated from the Chinese originals by an AI agent (Claude Code + glm-5.3); all formulas and code are kept verbatim from the Chinese version. The Chinese documents remain the authoritative version.
:::

:::{note}
Unless stated otherwise, each document treats only the closed-shell case; open-shell cases require the corresponding expansion over spin components (see the introductions of each document). Tensor shapes follow the column-major convention; the `becke-deriv` document follows the row-major convention in line with its implementation style, as explained at its beginning.
:::

```{toctree}
:maxdepth: 1
def
skeleton2
vmat1
becke-deriv
becke-grid-shift
```
