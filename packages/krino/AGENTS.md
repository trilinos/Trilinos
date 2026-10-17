# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## What Krino Is

Krino is the SIERRA package for **level-set geometry and conforming mesh
decomposition**. Its core job is **CDFEM** (Conformal Decomposition Finite
Element Method): given one or more implicit surfaces (level sets, analytic
surfaces, or faceted geometry), it cuts a background mesh so that element faces
align with the interfaces, producing phase-conforming subelements that physics
packages can solve on. It also provides the supporting machinery around that:

- **Level sets** — signed-distance fields, redistancing (`Redistance_Method`:
  `CLOSEST_POINT`, `FAST_MARCHING`, `FAST_ITERATIVE`), semi-Lagrangian advection.
- **Conformal decomposition** — element cutting, interface geometry, snapping,
  phase support, prolongation of solution data onto the new mesh.
- **Adaptive refinement** — h-adaptivity driven by the interface, plus general
  mesh refinement and refinement-aware rebalancing.
- **Mesh quality & smoothing** — quality metrics, mesh smoothing/optimization
  (some via Trilinos ROL).
- **Mesh utilities** — bounding-box mesh generation, mesh cloning, surface/facet
  handling.

Krino runs standalone (the `krino` executable, driven by a YAML input deck) and
is also consumed by other SIERRA apps — notably **Aria**, which links
`krino::krino_lib`, `krino::krino_sierra`, and `krino::krino_rebalance_utils` and
consumes the `krino_commands.xmldb` command DB to drive CDFEM (see
@../aria/AGENTS.md). Aria integrates through the `krino_sierra` framework/parser
bridge, not by including `Akri_*` headers directly.

## Building and Testing

See the top-level @../AGENTS.md for the `build_and_test` MCP tool rules.
Krino-specific selectors:

- Build only: `packages: ["krino"]` with `test_mode: "skip"`.
- Unit tests: building `krino` produces a **per-library** `krino_<lib>_utest`
  GTest executable (e.g. `krino_interface_geometry_utest`, `krino_decomp_utest`)
  rather than one monolithic `krino_unit`; use `test_mode: "gtest"` with
  `gtest_filter`. The unit binaries run under MPI — set `gtest_num_procs` to 1–4
  for distributed paths.
- Regression tests: `test_mode: "selector"` with `test_glob_pattern: "krino"`.

## Build Layout & CMake

**Two nested CMake levels.** The package root `krino/CMakeLists.txt` builds
either as a Trilinos/TriBITS package (`TRIBITS_PACKAGE(Krino)`) or as a
standalone CMake project (`SierraCmakeCode.cmake`). It adds two subdirectories:
`krino/` (the real code) and `delete_small_elements/` (a small standalone tool
exe, `Akri_DeleteSmallElementsMain.cpp`).

`krino/krino/CMakeLists.txt` is the hub: it `add_subdirectory`s each component
library, then builds the `krino` executable from `Apps_krino.cpp`.

**Sources are globbed.** Every component `CMakeLists.txt` uses
`FILE(GLOB HEADERS *.hpp)` / `FILE(GLOB SOURCES *.cpp)`, so adding a `.hpp`/`.cpp`
to an existing component directory needs no CMake edit — but a brand-new
component directory must be added to `krino/krino/CMakeLists.txt` (and given its
own `TRIBITS_ADD_LIBRARY`).

**Path-less includes.** Krino includes are angle-bracket and directory-less
(`#include <Akri_Foo.hpp>`); every build system puts each component's source dir
on the include path. So **moving a `.hpp`/`.cpp` between component dirs needs no
`#include` edits** in the file or its consumers — only the build-system library
membership and dependency edges change (plus the file's own include guard, which
encodes its path).

### Three build systems to maintain for krino

Although the top-level @../AGENTS.md says to ignore `bjam`, krino's own
`Jamfile` **is** kept in sync in practice. A change to krino's *library
structure* — adding, removing, moving, or re-wiring a component library — must be
mirrored in **all three** build systems, or one of the builds breaks:

| System | File(s) | Add-a-library primitive |
|---|---|---|
| **TriBITS CMake** (Trilinos build) | per-dir `krino/krino/<dir>/CMakeLists.txt` + `add_subdirectory` in `krino/krino/CMakeLists.txt` | `TRIBITS_ADD_LIBRARY(... DEPLIBS ...)` |
| **SierraCMake** (standalone) | `krino/SierraCmakeCode.cmake` | `krino_add_library(name krino/<dir> PUBLIC_<dep> ...)` |
| **bjam** (legacy) | `krino/Jamfile` | a `lib <name>` block + add to consumers' deps and the install list |

Checklist for a **new** component library (what the `field_calcs` addition
touched):
1. TriBITS: create `<dir>/CMakeLists.txt` with its `TRIBITS_ADD_LIBRARY` +
   `INSTALL`, and add `add_subdirectory(<dir>)` in `krino/krino/CMakeLists.txt`
   **before** its consumers.
2. SierraCMake: add a `krino_add_library(...)` block in `SierraCmakeCode.cmake`
   (deps use the `PUBLIC_`/`PRIVATE_`/`INTERFACE_` prefix convention); the
   `krino_add_library` helper auto-installs and registers the target in the
   `krinoTargets` export — no separate install/export edit needed.
3. bjam: add a `lib <name>` block in `Jamfile` and list it in each consumer's
   deps and in the user-install lib list.
4. Add the new library to every **consuming** library's dep list in all three
   systems. A leaf lib pulled in transitively (e.g. through
   `interface_geometry`/`decomp`) need not be added to the `krino_lib` umbrella.

**`krino_lib` is external-only.** The `krino_lib` umbrella exists solely for
out-of-repo consumers (Aria/platoengine) and the in-repo `pykrino` Python
binding. No target **inside** krino may depend on `krino_lib` in any of the three
build systems — internal consumers depend on the specific component libraries
they actually use (relying on `decomp`'s transitive closure where applicable).
Depending on the umbrella internally hides real dependencies and links every
target to everything, undermining the layering. When a build fails with a
missing symbol/header after removing an umbrella dep, add the specific owning
component lib — never re-add `krino_lib`.

Verify with the `build_and_test` MCP tool against the SierraCMake build
(`packages: ["krino"]`).

### Library targets (one per component dir under `krino/krino/`)

| Target | Dir | Contents |
|---|---|---|
| `krino_config_lib` | `config/` | build/config headers |
| `krino_diagwriter_lib` | `diagwriter/` | diagnostic/trace logging |
| `krino_master_element_lib` | `master_element/` | master (reference) element / shape functions |
| `krino_geometry_lib` | `geometry/` | geometric primitives, transforms |
| `krino_math_utils_lib` | `math_utils/` | math/linear-algebra helpers |
| `krino_quality_metric_lib` | `quality_metric/` | element quality metrics |
| `krino_quality_metric_sens` | `quality_metric_sens/` | quality-metric sensitivities |
| `krino_mesh_quality_lib` | `mesh_quality/` | mesh-wide quality evaluation |
| `krino_rol_lib` | `rol/` | Trilinos ROL optimization glue |
| `krino_smoothing_lib` | `smoothing/` | mesh smoothing / node optimization |
| `krino_surface_lib` | `surface/` | analytic/implicit surfaces |
| `krino_mesh_surface_lib` | `mesh_surface/` | faceted / mesh-based surfaces, incl. building interface-conforming facets from mesh sides (`FacetsFromSides`, `OrientedSideNodes`); `DEPLIBS krino_surface_lib krino_mesh_utils_lib` |
| `krino_mesh_utils_lib` | `mesh_utils/` | mesh construction/cloning utilities |
| `krino_refinement_lib` | `refinement/` | adaptive mesh refinement |
| `krino_refinement_rebalance_lib` | `refinement_rebalance/` | refinement-aware load rebalancing |
| `krino_field_calcs_lib` | `field_calcs/` | nodal field-computation utilities (bounding box, gradient) |
| `krino_support_lib` | `support/` | CDFEM/refinement config substrate (`CDFEM_Support`, `RefinementSupport`, `CDFEM_Snapper`, `LS_Field`) — below both `levelset` and `interface_geometry` |
| `krino_levelset_lib` | `levelset/` | level-set core: `LevelSet`, redistancing, semi-Lagrangian advection, surface distance, IC algorithms, surface registry (`DEPLIBS krino_support_lib` + contour/eikonal/surface/…) |
| `krino_interface_geometry_lib` | `interface_geometry/` | CDFEM cutting machinery: element cutters, parent edges, intersection points, `*InterfaceGeometry` classes (`DEPLIBS krino_levelset_lib krino_support_lib …`) |
| `krino_decomp_lib` | `decomp/` | top-level decomposition driving (`DEPLIBS krino_interface_geometry_lib …`) |
| `krino_lib` | `krino_lib/` | source-free **umbrella** re-exporting decomp, interface_geometry, post_process, mesh_clone, levelset_shape_sens, etc. **External-only** — see the external-only rule below; internal targets must depend on specific component libs, never this umbrella. |
| `krino_region_lib` | `region/` | Region/Simulation orchestration (`DEPLIBS krino_decomp_lib krino_post_process_lib`) |
| `krino_rebalance_utils_lib` | `rebalance_utils/` | load-balance utilities (`DEPLIBS krino_decomp_lib`) |
| `krino_parser_lib` | `parser/` | input-deck parsing (`DEPLIBS krino_region_lib krino_decomp_lib krino_post_process_lib krino_surface_lib krino_mesh_surface_lib`) |

### Executable / test targets
- `krino` — the production executable (`Apps_krino.cpp`).
- **Per-library GTest binaries** — each component library with a
  `unit_tests/` subdir builds its own `krino_<lib>_utest` executable (e.g.
  `krino_interface_geometry_utest`, `krino_decomp_utest`, `krino_surface_utest`),
  each linking `krino_unit_test_utils` (which pulls in the `krino_unit_main`
  gtest-main from `unit_main/`). There is no longer a single monolithic
  `krino_unit` binary. A library's tests link that library plus its transitive
  deps, so tests for relocated code run under whichever `_utest` owns their
  `unit_tests/` dir. The level-set-owned tests (`Akri_Unit_SemiLagrangian`,
  `Akri_Unit_OutputUtils`, `Akri_Unit_Constrained_Redistance`) live in
  `levelset/unit_tests/` and run in `krino_levelset_utest`. Note that
  `Constrained_Redistance` also drives `refinement` as test scaffolding;
  `refinement` and `levelset` are independent sibling libraries (neither depends
  on the other), and the `_utest` links `refinement` because
  `krino_unit_test_utils` now depends on `krino_refinement` **directly** (not via
  the `krino_lib` umbrella), so this is fine.
- Every per-library `krino_<lib>_utest` must be registered as a
  `TRIBITS_ADD_TEST_DIRECTORIES(<lib>/unit_tests)` in `krino/CMakeLists.txt`
  (after `unit_test_utils`) **and** have a matching
  `<lib>_unit_tests.test` in `tests/krino_rtest/krino/unit_tests/` (the runtime
  test suite) with `<Parameters exec="krino_<lib>_utest" .../>`. A test
  executable is a DAG leaf: it links the specific libraries it needs directly
  (e.g. `krino_math_utils_utest` links `krino_rol` for `Akri_ROLOptimize.hpp`,
  `krino_decomp_utest` links `krino_mesh_clone`), never the `krino_lib` umbrella.
- `delete_small_elements` — standalone mesh-cleanup tool.

## Framework Integration (`krino_sierra/`)

`krino_sierra/` (a sibling of `krino/` under the repo-root `krino/`) is the bridge
into the full SIERRA framework: `sierra::Sctl` events (`Akri_Events`,
`Akri_LevelSet_Sctl`), the `RegionInterface` singleton, `*_SierraParser` classes,
and the `xml/Akri_Levelset.xml` command grammar (built into `krino_commands.xmldb`).
This is what Aria links against to drive CDFEM. The related `krino_mesh_adapt/`
sibling is a standalone mesh-adaptation executable (`KrinoMeshAdaptMain.cpp` +
`mesh_adapt_lib/`, `Krino*` naming) that calls `krino::refine_mesh_with_params(...)`.

## Naming Conventions

- Source/header files use the **`Akri_`** prefix (e.g. `Akri_CDMesh.cpp`,
  `Akri_LevelSet.hpp`); headers are `.hpp`, sources `.cpp`. Exception:
  `krino_mesh_adapt/` uses `Krino*` names (`KrinoMeshAdapt*`).
- Everything lives in **`namespace krino`**. Classes inside are generally not
  prefixed (`LevelSet`, `CDMesh`, `Simulation`, `Region`).
- Unit-test sources also use `Akri_` — typically `Akri_Unit_*.cpp` for tests and
  `Akri_*Fixture.hpp` / `Akri_StkMeshBuilder` / `Akri_MeshSpecs.hpp` for
  fixtures/helpers.
- Toolkit siblings you'll see: `stk::mesh::{BulkData,MetaData}` (meshing
  infrastructure) and `sierra::Sctl` (framework, in `krino_sierra/` only).

## Code Style

Match the style of the file and directory you are editing — the conventions
below describe the prevailing patterns, but local consistency wins over any
global rule.

### Formatting
- **2-space indentation**, spaces not tabs.
- **Braces on their own line** (Allman): the opening brace of a function,
  class, `if`/`for`/`while`, etc. goes on the next line, not trailing the
  statement.
- Single-line `if`/loop bodies are commonly written without braces.

### Includes
- `#include` directives use angle brackets, not quotes (e.g.
  `#include <Akri_CDMesh.hpp>`), including for headers in the same directory.
- Group as: this file's own header first, then `<array>`/`<vector>`-style
  standard headers, then krino `Akri_*` and toolkit (`stk_*`, `Teuchos_*`,
  `ROL_*`) headers.

### Naming
- **Functions and methods**: `snake_case` (`compute_exact_winding_number`,
  `get_centroid`).
- **Local variables and parameters**: `camelCase` (`queryLoc`, `facetIndices`).
- **Types** (classes, structs, enums, aliases): `PascalCase`
  (`WindingNumberBottomUpTree`, `FacetClusterApproximation`).
- **Member variables**: `my`-prefixed camelCase (`myFacets`, `myRegions`,
  `myOps`). Existing code also contains older `my_snake` (`my_regions`) and
  `m_` (`m_facets`) forms; prefer `myVar` for new members, and match the
  prevailing prefix when editing a file that consistently uses another.
- Prefer well-named identifiers over comments (see below).

### Comments
- Prefer expressing intent through names over inline comments for routine code.
- Reserve comments for non-obvious rationale: policy/CRTP contracts, lifetime
  and ownership constraints, and tuning-parameter justifications (see the
  header comments in `geometry/Akri_WindingNumberFacetOps.hpp` for the expected
  level of a contract comment).

### Files
- BSD-style Sandia/NTESS license header at the top of every source file.
- Headers use `#ifndef KRINO_..._HPP_` include guards (not `#pragma once`).
- See Naming Conventions above for the `Akri_` file-prefix rule.

## Architecture

### Control flow
`Apps_krino.cpp` `main()`:
1. `krino::Startup` (`region/Akri_Startup`) — boots MPI / the runtime
   environment; `exit_early()` handles `--help`-style early exits.
2. `krino::Simulation simulation(...)` (`region/Akri_Simulation`) — the
   top-level owner (`std::vector<std::unique_ptr<Region>> my_regions`).
3. `krino::Parser::parse(simulation)` (`parser/Akri_Parser`, sub-namespace
   `krino::Parser`) — reads the **YAML** input deck (guarded by
   `KRINO_HAVE_YAML`) and populates the simulation with regions, level sets,
   phases, CDFEM options, etc. via per-feature parsers in `parser/`.
4. `simulation.commit()` → `initialize()` (currently an inline no-op) →
   `execute()` — the standard SIERRA commit/initialize/run lifecycle. Exceptions
   are funneled through `Startup::handle_exception`.

### The Region / Simulation (`region/`)
`Simulation` owns one or more `Region`s; a `Region` owns a mesh, its fields, its
level sets, and the decomposition/adaptivity algorithms, and advances the
solution. This is the central orchestrator — most cross-cutting behavior lands
here.

### Core decomposition engine (`krino_lib/`) — the largest, most important dir
- **`Akri_CDMesh`** — the conforming decomposition mesh: the heart of CDFEM.
  Manages the background mesh, its decomposition into phase-conforming
  subelements, and prolongation of data across a re-decomposition.
- **Interface geometry** — `Akri_InterfaceGeometry` (abstract),
  `Akri_LevelSetInterfaceGeometry`, `Akri_AnalyticSurfaceInterfaceGeometry`,
  `Akri_CreateInterfaceGeometry`. Describes the surfaces to cut against.
- **Cutting** — `Akri_Element_Cutter`, `Akri_Cutting_Surface`,
  `Akri_SubElement`, `Akri_ContourSubElement`, `Akri_ContourElement`,
  `Akri_AdaptiveContourTet/Tri`, `Akri_CDFEM_Parent_Edge(s)`.
- **Phases** — `Akri_Phase_Support`, `Akri_ConformingPhaseParts`,
  `Akri_CDFEM_Support`, `Akri_CDFEM_Snapper` — map interfaces/level-set signs to
  mesh parts (phases) and control snapping tolerances.
- **Level sets** — `Akri_LevelSet` (one per captured interface), redistancing
  via `Akri_ClosestPointRedistance`, `Akri_Fast_Marching`,
  `Akri_FastIterativeMethod` (`Akri_Eikonal_Calc`, `Akri_DistanceSweeper`);
  `Akri_SemiLagrangian` advection; `Akri_LevelSetShapeSensitivities` for
  topology-optimization sensitivities.
- **Snapping & prolongation** — `Akri_Snap`, `Akri_ProlongationData`,
  `Akri_MeshClone`, `Akri_DecompositionHasChanged`.
- **Refinement glue** — `Akri_RefinementManager`, `Akri_CDMesh_Refinement`.
- **Mesh access** — `Akri_AuxMetaData` wraps STK meta/field registration.

### Refinement (`refinement/`, `refinement_rebalance/`)
Adaptive h-refinement of the STK mesh, driven by the interface or error
indicators, with a parallel-rebalance step that preserves parent/child
refinement relationships.

### Surfaces, quality, smoothing
- `surface/` + `mesh_surface/` — implicit/analytic surfaces and faceted
  (mesh-based) surfaces used as CDFEM interface geometry.
- `quality_metric/` + `quality_metric_sens/` + `mesh_quality/` — element/mesh
  quality metrics and their sensitivities.
- `smoothing/` + `rol/` — mesh smoothing / node-position optimization, some
  formulated as ROL optimization problems.

## Unit Test Patterns

Tests are GTest sources in each library's `<lib>/unit_tests/Akri_Unit_*.cpp`,
built into that library's `krino_<lib>_utest` binary. Fixtures/helpers build
small meshes analytically:
`Akri_StkMeshFixture.hpp`, `Akri_StkMeshBuilder`, `Akri_MeshFromFileFixture.hpp`,
`Akri_MeshSpecs.hpp` / `Akri_SecondOrderMeshSpecs.hpp`. Tests construct a mesh,
apply a level set / analytic surface, run a cutter / decomposition / refinement /
smoothing step, and check the resulting subelements, phases, or quality. Heavier
fuzz/random/performance cases are named so they can be included or excluded via
`gtest_filter` (e.g. `*fuzz*`, `*Random*`, `*Refinement.performance*`).

### Pitfall: `-Werror=dangling-else` with streamed gtest assertions

A *streamed* gtest assertion — `EXPECT_*`/`ASSERT_*` with a trailing
`<< "message"` — expands to an internal `if (...) ... else ...`. When such an
assertion is the unbraced body of your own `if`, the compiler sees a dangling
`else` and the CI build fails with `-Werror=dangling-else`. Always wrap the
body in braces in this case, even though single-line `if` bodies are otherwise
written without braces here:

```cpp
if (mesh.parallel_size() > 1)
{
  EXPECT_TRUE(cond) << "message";   // braces required: streamed assertion under an if
}
```

An unbraced body is fine for a non-streamed assertion (no `<<`).

Only gcc reports this; the default `docker-clang` build of the `build_and_test` tool does not, so
a clean clang build does not show the code is free of it.
