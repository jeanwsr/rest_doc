# Constants Module (constants)

## Overview

The constants module centralizes the global constants used throughout the REST codebase, including physical and mathematical constants, numerical thresholds, the index constants required by the libcint interface, and periodic-table-related data. The module entry point is `src/constants/mod.rs`, with several submodules holding larger data tables.

Centralizing constants serves to:

- **Single source of truth**: each physical quantity has exactly one definition across the codebase, avoiding inconsistent hardcoded duplicates;
- **Traceability**: physical constants carry source annotations (e.g. the CODATA year) for easy verification and updating;
- **Ease of maintenance**: a single update propagates to all usage sites.

## Module Structure

| Submodule | File | Contents |
|---|---|---|
| `mod` (root) | `src/constants/mod.rs` | Physical and mathematical constants, unit conversion constants, numerical thresholds, scale constants, libcint index constants |
| `element` | `src/constants/element.rs` | Element names, atomic masses and charges, electron configurations, atomic radii, shell information |
| `vsap` | `src/constants/vsap.rs` | Numerical potential tables for the VSAP (superposition of atomic potentials) initial guess |
| `c2s` | `src/constants/c2s.rs` | Cartesian-to-spherical basis transformation matrices |
| `cartesian_gto` | `src/constants/cartesian_gto.rs` | Layout information for Cartesian Gaussian-type orbitals (GTOs) |
| `solvent` | `src/constants/solvent.rs` | Solvent parameter data used by solvent models |

Each submodule is re-exported in the root module via `pub use`, so all constants can be referenced uniformly from `crate::constants`.

## Physical and Mathematical Constants

The physical constants are defined in `src/constants/mod.rs`, each with a source annotation. Note that these constants originate from **different CODATA release years** and are not fully internally consistent (see the "Source" column below).

| Constant | Value | Meaning | Source |
|---|---|---|---|
| `EV` | 27.2113845 | Hartree → eV | CODATA 2002 |
| `FQ` | 1822.8884861920776 | Ratio of atomic mass unit u to electron mass mₑ (= 1 / electron mass in u) | CODATA 2014 |
| `E` | `std::f64::consts::E` | Euler's number e | mathematical constant |
| `PI` | `std::f64::consts::PI` | π | mathematical constant |
| `LIGHT_SPEED` | 137.03599967994 | Inverse fine-structure constant α⁻¹ (speed of light in atomic units) | CODATA 2006 |
| `BOHR` | 0.52917721092 | Bohr radius (in Å), used for Å↔Bohr conversion | CODATA 2010 |
| `BOHR_SI` | `BOHR * 1e-10` | Bohr radius (in m) | derived |
| `G_ELECTRON` | 2.00231930436182 | Electron g-factor | CODATA 2014 |
| `E_MASS` | 9.10938356e-31 | Electron mass (kg) | CODATA 2014 |
| `AVOGADRO` | 6.022140857e23 | Avogadro constant | CODATA 2014 |
| `PLANCK` | 6.626070040e-34 | Planck constant (J·s) | CODATA 2014 |
| `BOLTZMANN` | 1.380649e-23 | Boltzmann constant (J/K) | CODATA 2018 (exact, SI 2019) |
| `CLIGHT_CMS` | 2.99792458e10 | Speed of light (cm/s) | SI definition (exact) |
| `R_GAS` | `BOLTZMANN * AVOGADRO` | Ideal gas constant (J/(mol·K)) | derived |
| `E_CHARGE` | 1.6021766208e-19 | Elementary charge (C) | CODATA 2014 |
| `DEBYE` | 3.335641e-30 | 1 debye in C·m | defined via the speed of light |

```{note}
For historical reasons, the constants above are taken from several different CODATA releases, most of them (`G_ELECTRON`, `E_MASS`, `AVOGADRO`, `PLANCK`, `E_CHARGE`) from CODATA 2014. Newer constants (`HARTREE2J`, `HARTREE2WAVENUMBERS`, `AMU2KG`, `SPEED_OF_LIGHT`) use CODATA 2022. If you wish to unify them to a single release, please carefully assess the impact on regression test baselines.
```

## Unit Conversion Constants

These constants handle conversions between different unit systems. They are also defined in `src/constants/mod.rs`, grouped under the `unit conversion` comment block.

| Constant | Value | Meaning | Source |
|---|---|---|---|
| `HARTREE2J` | 4.359744722206e-18 | Hartree → joule | CODATA 2022 |
| `HARTREE2KJMOL` | `HARTREE2J * AVOGADRO / 1000.0` | Hartree → kJ/mol (≈2625.5) | derived |
| `HARTREE2KCALMOL` | `HARTREE2KJMOL / CALORIE2J` | Hartree → kcal/mol | derived |
| `HARTREE2KCAL` | 627.5094841703362 | Hartree → kcal/mol (= HARTREE2KCALMOL) | derived |
| `HARTREE2WAVENUMBER` | 219474.63 | Hartree → cm⁻¹ | CODATA |
| `HARTREE2WAVENUMBERS` | 219474.63136314 | Hartree → cm⁻¹ (higher precision) | CODATA 2022 |
| `CALORIE2J` | 4.184 | calorie → joule | Wikipedia |
| `AMU2KG` | 1.66053906892e-27 | amu → kg | CODATA 2022 |
| `SPEED_OF_LIGHT` | 299792458.0 | Speed of light (m/s) | CODATA 2022 |
| `BOHR2ANG` | `BOHR` | Bohr → Å (alias) | alias |
| `AU2DEBYE` | `E_CHARGE * BOHR*1e-10 / DEBYE` | Atomic-unit dipole → debye (≈2.541746) | derived |

## Numerical Thresholds and Scale Constants

| Constant | Value | Meaning |
|---|---|---|
| `AUXBAS_THRESHOLD` | 1.0e-10 | Threshold for the inverse (and inverse square root) of the auxiliary Coulomb matrix |
| `INVERSE_THRESHOLD` | 1.0e-10 | Matrix inversion threshold |
| `SQRT_THRESHOLD` | 1.0e-10 | Matrix square-root threshold |
| `E5` – `E9` | 1.0e5 – 1.0e9 | Common order-of-magnitude scale factors |
| `MPI_CHUNK` | 134217728 | MPI transfer chunk size (about 1 GB) |

## libcint Index Constants

These integer constants are used to access fixed slots in libcint's `bas`, `atm`, and `env` arrays, and must be kept consistent with libcint's conventions.

| Category | Constants |
|---|---|
| `bas` array indices | `BAS_ATM`, `BAS_ANG`, `BAS_PRM`, `BAS_CTR`, `BAS_SLOTS` |
| `atm` array indices | `ATM_NUC`, `ATM_ENV`, `ATM_NUC_MOD_OF`, `ATM_FRAC_CHARGE_OF`, `ATM_SLOTS` |
| ECP-related | `ECP_LMAX`, `NUC_ECP` |
| `env` pointers and nuclear models | `PTR_EXPCUTOFF`, `PTR_COMMON_ORG`, `PTR_RINV_ORIG`, `NUC_MOD_OF`, `NUC_STAD_CHARGE`, `NUC_GAUS_CHARGE`, `NUC_FRAC_CHARGE`, `ENV_PRT_START` |

For more details on libcint data structures, refer to the code agent skill `rest-libcint-knowledge` (if available).

## Element Data (`element`)

Periodic-table-related data is collected in `src/constants/element.rs`:

| Constant / static | Type | Meaning |
|---|---|---|
| `SPECIES_NAME` | `[&str; 118]` | Element symbols (H to Og) |
| `MASS_CHARGE` | `[(f64, f64); 118]` | (atomic mass, nuclear charge) for each element; atomic masses use the IUPAC 2021 standard atomic weights |
| `ELEM1ST` – `ELEM7TH` | `[&str; N]` | Element symbol lists grouped by period |
| `ELEMTMS` | `[&str; 40]` | List of transition-metal elements |
| `ATOM_CONFIGURATION` | `[[usize; 4]; 119]` | Electron filling configuration of each element in the s/p/d/f shells |
| `S_SHELL`/`P_SHELL`/`D_SHELL`/`F_SHELL` | `[f64; N]` | Occupation templates for each angular-momentum shell (used by SAD and ECP) |
| `XE_SHELL`/`KR_SHELL` | `[f64; N]` | Shell occupation templates for noble-gas cores |
| `NELE_IN_SHELLS` | `[f64; 15]` | Number of electrons each shell can hold |
| `SPECIES_INFO` | `HashMap<&str, &(f64,f64)>` | Map from element symbol to (mass, charge) (`lazy_static`) |
| `ATOMIC_RADII` | `HashMap<&str, &f64>` | Map from element symbol to atomic radius (Å) (`lazy_static`) |

## Usage Conventions

To keep the codebase consistent, please follow these conventions:

1. **Import constants with `use` at the top of the file**; do not use fully qualified paths inside function bodies.

    ```rust
    // Recommended
    use crate::constants::BOHR;
    // ...
    let dist_bohr = dist_ang / BOHR;

    // Not recommended
    let dist_bohr = dist_ang / crate::constants::BOHR;
    ```

2. **Do not hardcode physical constant values** in the code (such as `0.52917721092`, `27.2114`, `1822.8885`, etc.). Always reference them from the `constants` module.

3. If a required constant does not yet exist, add it to `src/constants/mod.rs` with a source annotation (e.g. the CODATA year and the corresponding NIST value).

## Sources and Updates

The values of the physical constants are taken from the [NIST CODATA](https://physics.nist.gov/cuu/Constants/) database. Historical releases can be verified against the archive files, for example:

```
https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2010.txt
```

When updating a constant, indicate the adopted CODATA release year in its inline comment, and record the value change in the commit message so that its impact on computed results (energies, gradients, frequencies, etc.) and regression baselines can be assessed.
