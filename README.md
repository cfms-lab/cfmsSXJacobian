# cfmsSXJacobian — automatically derived bending forces and Jacobians for implicit garment simulation

This repository accompanies the manuscript

> InHwan Sul, *Effects of approximate bending derivatives in implicit garment simulation: a spurious edge force and a numerical drag*,
> submitted to The Visual Computer (2026).

It contains a code generator that takes one condition function `C(x)` of a Baraff–Witkin style force, derives `C`, `∇C` and `∇²C`
symbolically (SymPy) and writes C++ evaluation and assembly functions; the specifications of the forces used in the manuscript; the
tests that verify the generated derivatives against Mathematica and finite differences; the measurement harnesses behind the tables
and figures; and the recorded results of every run.

## About this snapshot

- This is a **single-commit snapshot** of the private development repository (commit `8b37066`, 2026-11-15). The history is not included.
- Not included: the manuscript sources and figure scripts (`draft/`), the third-party garment files of the collar measurement
  (`data/garments/`, see [Garment data](#garment-data)) and internal working notes.
- The design and judgment documents in `docs/` and most comments are written in **Korean**; `README.ko.md` is the Korean overview.
  Each document states its acceptance criteria and predictions, written before the run, and the outcome.
- **Units.** The result files use the units of the drape engine (CGS: cm, g, s; erg, dyn/cm, g/cm²). The manuscript uses SI:
  1 erg = 10⁻⁷ J, 1 dyn/cm = 10⁻³ N/m, 1 g/cm² = 10 kg/m², speeds in cm/s are divided by 100.

## Layout

| Path | Contents |
|---|---|
| `codegen/sxgen.py` | Generator: specification → `generated/<NAME>.h` (`Eval`, `EvalFrozen`, `EvalStaged`, `EvalImplicit`, `EvalImplicitGN`, …) and the engine-form header `generated/<NAME>_Engine.h` |
| `codegen/specs/*.py` | Specifications (condition function, node and stencil variables). Bending: `baraff_bend.py` (C = θ), `bridson_bend.py` (C = 2√(2w̄) sin(θ/4)) |
| `codegen/sxview.py`, `notebooks/` | Tensor-form view of the derivatives (Jupyter, prototype; independent of the generator) |
| `src/` | Force and Jacobian assembly from `C`, `∇C`, `∇²C` (`sx_condition_force.h`); variables and stencil gathering (`sx_stencil.h`) |
| `third_party/sxmesh/` | SXMesh (named structure-of-arrays container) — see [License](#license) |
| `model/` | Test samples and independent references (Mathematica `.wl`) |
| `tests/` | Gate tests (C++ and Python), measurement harnesses, diagnostics (`diag_*`) |
| `results/` | Recorded results (`*.json`) and Mathematica references (`mma_ref*.txt`) |
| `docs/` | Design and judgment documents (Korean) |
| `run_gate.bat` | Builds and runs the gates |

`generated/` is not tracked; it is rebuilt deterministically from the specifications.

## Requirements

- Windows (the gate script is a batch file), Python 3 with **SymPy 1.14** and **NumPy**, **MSVC 2022** (`vcvars64`; the code is C++11).
- Optional: Mathematica 14.3 (`wolframscript`) for the independent references. Without it, `run_gate.bat nomma` uses the committed
  references in `results/mma_ref*.txt`.
- The engine measurements need the CPU build of the cfmsDrape drape engine (cfms-lab), which is **not public**. See [Engine measurements](#engine-measurements).

## Running

```
python codegen\sxgen.py codegen\specs\bridson_bend.py --out generated --results results
run_gate.bat nomma
```

`run_gate.bat` generates the C++ from the specifications, builds the tests (output in `%TEMP%\cfmsSXJacobian_build`) and runs the gates;
the numbers are written to `results\*.json`. Without the engine, the engine-dependent stages are skipped with a message, and the
specification, Mathematica, finite-difference, assembly, integration and timing gates run. Some gates were closed as NO-GO by their
pre-set criteria (for example the assembly gate, `docs/GATE_ASSEMBLY.md`); the summary line lists them. Timing values depend on the machine.

## Manuscript and repository

The figures were drawn from the result files below by scripts that are not part of this snapshot.

| Manuscript | Document | Harness | Result | Engine needed |
|---|---|---|---|---|
| Table 1 (verification of the generated derivatives) | `docs/MODEL_PROBLEM.md`, `docs/GATE_BRIDSON.md` | `tests/test_baraff_bend.cpp`, `tests/test_bridson.cpp` (via `run_gate.bat`) | `results/gate_baraff_bend.json`, `results/gate_bridson.json` | no |
| Sect. 4.2, Table 2, Fig. 3, Table S1 (spurious edge force on a corrugated sheet) | `docs/MEASURE_FROZEN_FORCE.md` (M1–M3) | `tests/measure_frozen_force.cpp` | `results/measure_frozen_force.json` | no |
| Fig. 4, Table S2 (collar, linear response) | `docs/MEASURE_FROZEN_FORCE.md` (M4b) | `tests/measure_frozen_garment.cpp`, `tests/measure_frozen_garment_clo.py` | `results/measure_frozen_garment_clo.json` | no — garment files needed |
| Table S3 (collar, nonlinear equilibrium) | `docs/MEASURE_FROZEN_FORCE.md` (M4c) | `tests/measure_frozen_garment_nl.cpp`, `.py` | `results/measure_frozen_garment_nl.json` | no — garment files needed |
| Table 3, Figs. 5 and 6, Tables S4 and S5 (one-mode step, terminal velocity) | `docs/GATE_JACOBIAN_DAMPING.md` | `tests/gate_jacobian_damping.py`, `_2.py`, `tests/measure_mode_sum.py`, `tests/diag_block_symmetric_kappa.py` (formula only) | `results/gate_jacobian_damping.json`, `_2.json`, `results/measure_mode_sum.json` | yes |
| Table 4, Figs. 7 and 8 (fall of a released sheet) | `docs/GATE_JACOBIAN_DAMPING.md` (section 「논문용 측정」) | `tests/measure_jacobian_drop.py` | `results/measure_jacobian_drop.json`, `results/fig_drop_render.json` | yes |
| Table S6 (six bending matrices) | `docs/GATE_JACOBIAN_DAMPING.md` (sections 「논문용 측정 3」, 「논문용 측정 4」) | `tests/engine_bend_exact.py`, `tests/measure_jacobian_drop_exact.py`, `tests/engine_bend_symgn.py`, `tests/measure_jacobian_symgn.py`, `tests/diag_symgn_hinge.py`, `tests/diag_symgn_folded.py` | `results/gate_engine_bend_exact.json`, `results/measure_jacobian_drop_exact.json`, `results/gate_engine_bend_symgn.json`, `results/measure_jacobian_symgn.json`, `results/diag_symgn_*.json` | yes (experimental engine modes) |
| Table S7 (static drape coefficient) | `docs/GATE_ENGINE_SUITE.md` (「판정 뒤 진단」, item 6) | `tests/diag_cusick_dc.py` | `results/diag_cusick_dc.json` | yes |
| Table S8 (small-step reference solution) | `docs/GATE_LINK_L4R.md`, `docs/GATE_LINK_L4.md` | `tests/link_l4r.py`, `tests/link_l4.py` | `results/gate_link_l4r.json`, `results/gate_link_l4_2.json` | yes |
| Table S9 (computation time, stiff bending) | `docs/GATE_JACOBIAN_DAMPING.md` (section 「논문용 측정 2」) | `tests/measure_time_stiff.py` | `results/measure_time_stiff.json` | yes |
| Sect. S7 (time per hinge) | `docs/GATE_GAUSS_NEWTON.md`, `docs/GATE_BRIDSON.md` | `tests/timing_bases.py` | `results/gate_timing_gn.json`, `results/gate_timing_bridson.json` | yes |

## Engine measurements

The harnesses marked "yes" above load the CPU DLL of the cfmsDrape engine through its Python wrapper and debug hooks. The engine is
not public; its build script is `tests/build_engine_dll.bat` (`CFMS_DRAPE_DIR` points to an engine checkout). The engine-side
changes that the measurements use — the generated Bridson Gauss–Newton bending kernel and the generated stretch and shear kernels,
and the experimental bending modes `exact_state`, `psd_state` and `symgn_state` — are described in `docs/GATE_LINK_L2.md`,
`docs/GATE_LINK_L3.md` and `docs/GATE_JACOBIAN_DAMPING.md`. Every run is recorded in `results/`, including the engine commit and
build state.

## Garment data

The collar measurements (Fig. 4, Tables S2 and S3) use three free items of **CLO-SET CONNECT** (creator: CONNECT Official). Their
licence does not allow redistribution, so the files are not included. They can be downloaded from the CLO-SET CONNECT website
(an account is needed):

| Folder in `data/garments/` | Item | Page | CLO version measured | SHA-256 of the `.zprj` measured |
|---|---|---|---|---|
| `kofoti_024_classic_dress_shirt` | MV2 Classic Dress Shirt | <https://connect.clo-set.com/detail/6ffbd199f8fc4ca48159082cbb5cc1d6> | 2025.0.128 | `7c3c2cc163a4318116b25ab3d343af9ef0291a596046ce23e2aec5973e9eca62` |
| `kofoti_370_button_down_shirt` | MV2 Button-Down Shirt | <https://connect.clo-set.com/detail/7a448111957c4967acb3807241c404d8> | 2024.2.120 | `bcb179cfbeb51f93cb6d341334068a86a2f9511fa41f4bd31ef962ea4dbc89ba` |
| `kofoti_244_wide_fit_open_collar_shirt` | MV2 Wide Fit Open Collar Shirt | <https://connect.clo-set.com/detail/81b1095d5db34164a73886e6f8f94b9d> | 2025.0.432 | `7621cbf4b407e469c56dbd1e07c30f004ef0e9c10e09b9d886911db1c988a6cd` |

An item may have been updated since the measurement; the SHA-256 shows whether a downloaded file is the one measured.

1. Put each `.zprj` file in `data/garments/<folder>/source/`.
2. `python tests\export_clo_garment.py` unpacks each project into `garment.npz` (2D pattern mesh, the 3D positions draped by CLO,
   seams) and `garment.json`, and writes `materials.json`. This step uses a reader for the CLO project format that is not part of
   this repository (option `--testbed`).
3. Build `tests/measure_frozen_garment.cpp` (`docs/MEASURE_FROZEN_FORCE.md`, section 「재현」) and run
   `python tests\measure_frozen_garment_clo.py` (Fig. 4, Table S2). Table S3: `tests/measure_frozen_garment_nl.cpp` and `.py`.

Collar patterns and properties used (`materials.json`; the other patterns are not simulated):

| Item | Collar pattern | Stretch `E` | Shear `SE` | Bending `BE` |
|---|---|---|---|---|
| all three | `10`, `8`, `6` (in the order of the table above) | 10⁶ dyn/cm (10³ N/m) | 5 × 10⁴ dyn/cm (50 N/m) | 10, 100, 1,000, 10,000 dyn·cm (1, 10, 100, 1,000 µJ) |

The numbers derived from the garments are recorded in `results/measure_frozen_garment_clo.json` and
`results/measure_frozen_garment_nl.json`.

## License

MIT (`LICENSE`), with the scope note in that file: `third_party/sxmesh/` and `tests/engine_ref/SDEPreHeaders.h` come from the
cfmsDrape engine (cfms-lab) and are not covered by the MIT licence. The garment items are subject to the CLO-SET CONNECT terms.

## Contact

InHwan Sul, Department of Materials Design Engineering, Kumoh National Institute of Technology, Gumi, Republic of Korea.
