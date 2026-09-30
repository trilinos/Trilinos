# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## What STK Is

STK (the **S**andia **T**ool**k**it) is the low-level infrastructure layer the
SIERRA apps are built on — most importantly the **parallel unstructured mesh**
(`stk_mesh`) plus the meshing/parallel utilities around it. It is a collection of
loosely-coupled modules (sub-packages); an app pulls in only the ones it needs
(e.g. `stk_search` + `stk_transfer` without `stk_mesh`). `stk_util` sits at the
bottom and depends on no other STK module.

**This repo is STK's home.** STK is developed here in the SIERRA tree and
*snapshotted into Trilinos periodically*, not the reverse. That dual-homing is why
STK looks different from the SIERRA-native packages: it uses the Trilinos
conventions — `.cpp`/`.hpp` extensions, the `stk::` namespace, and a
**TriBITS**-style sub-package/CMake system — and can build either standalone or as
a Trilinos package (`HAVE_STK_Trilinos`). When working here, treat SIERRA as the
context of record.

See the top-level @../AGENTS.md for the repo overview.

## Building and Testing

See the top-level @../AGENTS.md for the `build_and_test` MCP tool rules. STK-specific notes:

- Build only: `packages: ["stk"]`. STK builds many
  fine-grained library targets (`stk::stk_mesh_base`, `stk::stk_io`,
  `stk::stk_util_parallel`, `stk::stk_search`, ...) — downstream packages link
  the specific ones they use, not a single monolithic `stk` lib.
- Unit tests: each module has its own GTest executable named
  **`stk_<module>_utest`** (e.g. `stk_mesh_utest`, `stk_io_utest`,
  `stk_search_utest`), produced via TriBITS `TRIBITS_ADD_EXECUTABLE_AND_TEST`.
  Use `test_mode: "gtest"` with the matching executable/filter.

## Module Layout & Dependencies

Each module lives in a `stk_<name>/` directory with a doubled inner path
(`stk_mesh/stk_mesh/base/...`), its own `CMakeLists.txt` using `STK_SUBPACKAGE(...)`,
and a `cmake/Dependencies.cmake` declaring inter-module and TPL dependencies
(`LIB_REQUIRED_DEP_PACKAGES` / `..._DEP_TPLS`). The top-level `CMakeLists.txt`
drives everything through `STK_SUBPACKAGES()` / `STK_PACKAGE_POSTPROCESS()`
(macros in `cmake/stk_wrappers.cmake`).

The rough dependency stack (declared, not assumed — see each
`Dependencies.cmake`):
- **`stk_util`** — foundation: parallel communication (`CommSparse`, MPI
  wrappers), diagnostics/`DiagWriter`, command-line/env, reporting. Depends on no
  other STK module.
- **`stk_topology`** — element/side node-ordering and topology definitions
  (depends on `Shards`).
- **`stk_mesh`** — the parallel unstructured mesh: `BulkData` (entities +
  parallel connectivity + modification), `MetaData` (parts/fields/topologies),
  `Bucket`s (homogeneous entity storage), `Field`s. Requires `STKTopology`,
  `STKUtil`, `Kokkos`, and the `BLAS`/`MPI` TPLs.
- **`stk_io`** — reads/writes `stk_mesh` to Exodus via SEACAS
  (`SEACASIoss`/`SEACASExodus`); requires `STKMesh`.
- Peers used across the stack: **`stk_search`** (bounding-box proximity search,
  optionally accelerated by `ArborX`), **`stk_transfer`** (field transfer between
  meshes), **`stk_balance`** (Zoltan/parallel rebalance), **`stk_coupling`**
  (MPMD MPI coupling / `SplitComms`), **`stk_middle_mesh`** (common refinement of
  two surface meshes), **`stk_math`** / **`stk_simd`** (SIMD vector-instruction
  abstraction), **`stk_expreval`** (string expression evaluation), and the
  `*_util` companions (`stk_search_util`, `stk_transfer_util`, `stk_middle_mesh_util`).

## Naming & Conventions

- Sources/headers are `.cpp`/`.hpp` (Trilinos style) — **not** the `.C`/`.h` used
  by SIERRA-native packages.
- Everything is under `namespace stk`, with one sub-namespace per module:
  `stk::mesh`, `stk::io`, `stk::search`, `stk::topology`, `stk::balance`,
  `stk::coupling`, `stk::transfer`, `stk::util`, `stk::simd`,
  `stk::middle_mesh`, `stk::expreval`, plus `::impl` for internals.
- `STK_BUILT_FOR_SIERRA` gates SIERRA-specific capability; it also turns on
  `STK_ENABLE_16BIT_UPWARDCONN_INDEX_TYPE` (limits upward-connectivity buckets to
  ~65K entries — a deliberate SIERRA build choice, see the top-level
  `CMakeLists.txt`). In the SIERRA build this and the other STK CMake options
  (`STK_ENABLE_MPI`, `STK_ENABLE_ARBORX`, `STK_ENABLE_TESTS`, ...) are set by the
  Spack recipe's `cmake_args()`, not toggled by hand — see
  `spack_repo/sierra/packages/stk/package.py` in the environments repo (per the
  top-level @../AGENTS.md).

### NGP / performance portability
`stk_mesh` has a Kokkos device-mesh path for GPU execution: `DeviceMesh`,
`DeviceField`/`NgpField`, `DeviceBucket`, and the `GetNgpMesh` / `GetNgpField` /
`GetNgpExecutionSpace` accessors (the "NGP" = Next-Generation Platform naming used
across SIERRA). Kernels that must run on device go through these rather than
host-side `BulkData`/`Field` APIs. The device-mesh build is gated by the Spack
`+device_mesh` variant → `STK_USE_DEVICE_MESH` (with `+unified_memory` /
`+field_bounds_check` companions; `+simd` conflicts with `+cuda`/`+rocm`) — see
the stk `package.py`. `stk_ngp_test` is STK's GTest-like assertion framework that
works **inside device kernels** (`NgpTestDeviceMacros.hpp`), used by the NGP tests.

## Test Layout

Unlike the SIERRA-native packages (which keep tests in a per-package
`unit_tests/`), STK centralizes tests by module under top-level trees:
- **`stk_unit_tests/stk_<module>/`** — the GTest unit tests that build the
  `stk_<module>_utest` executables (`stk_unit_test_utils`, `stk_unit_main`,
  and `stk_mesh_fixtures` provide the shared harness/fixtures).
- **`stk_doc_tests/stk_<module>/`** — documentation/example tests: small,
  readable snippets that double as the code samples in the STK Manual.
- `stk_integration_tests/` and `stk_performance_tests/` — heavier end-to-end and
  benchmarking suites.

`CHANGELOG.md` tracks releases and `README.md` links the STK Manual PDF and lists
the modules authoritatively — consult them when in doubt about scope.
