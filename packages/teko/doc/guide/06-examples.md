# Examples {#teko_examples}

[TOC]

Teko ships four worked examples under `packages/teko/examples/`. This page summarizes what
each demonstrates and how to run it. The full sources are included in the API documentation
via the Doxygen `\include` directives at the bottom of each section.

> **Data files.** The drivers read Matrix-Market files from `examples/data/`, which ships
> `nslhs_test.mm` and `nsrhs_test.mm`. All three example programs (both `BuildPreconditioner`
> drivers and `StridedSolve`) additionally read `../data/nsjac.mm`, which is **not** currently
> in `examples/data/` — supply your own Navier–Stokes Jacobian there before running them.

All examples are built by the package's CMake when Teko's tests/examples are enabled; the
binaries land next to their source directories in the build tree.

---

## BuildPreconditioner — the canonical driver

`examples/BuildPreconditioner/` contains two drivers that read a strided Navier–Stokes
matrix, build a block preconditioner, and solve with Belos.

- **`example-driver.cpp`** — the *native Tpetra* path. It splits the monolithic matrix with
  a `StridedTpetraOperator` (`{2, 1}` = 2 velocity + 1 pressure per node), builds an LSC
  preconditioner via `InvLSCStrategy`, wraps it in a `TpetraBlockPreconditioner`, and hands
  it to `Belos::BlockGmresSolMgr` as a right preconditioner.
- **`example-driver-belos.cpp`** — the *Thyra/Stratimikos* path. It assembles a 2×2 block
  operator with `Thyra::block2x2`, builds an `InverseLibrary` from a hand-built parameter
  list (a `"Gauss-Seidel"` entry over Ifpack2 block solves), calls `Teko::buildInverse`, and
  solves with Belos. This is the driver dissected in
  [Getting Started](01-getting-started.md).

Run (from the example's build directory):

```console
$ mpirun -np 2 ./Teko_example-driver-belos.exe
```

**Full source** (Thyra/Stratimikos path):
\include example-driver-belos.cpp

**Full source** (native Tpetra path):
\include example-driver.cpp

---

## StridedSolve — SIMPLE on a strided operator, solver chosen at runtime

`examples/StridedSolve/strided_solve.cpp` reads `../data/nsjac.mm`, splits it with a strided
operator, and builds a **SIMPLE** preconditioner
(`Teko::NS::SIMPLEPreconditionerFactory(inverse, alpha)`). The inner solver name is read from
`solverparams.xml` (an ordinary Stratimikos solver list), so you can swap Amesos2/Belos/…
without recompiling. A good template for "drive the inner solver from XML while building the
block structure in code."

**Full source:**
\include strided_solve.cpp

---

## Arbitrary Tpetra blocking — field blocks from monolithic GIDs

The strided examples assume a regular nodal ordering such as `[u v p]` at every mesh node. For
application orderings that are not regular strides, use
`Teko::TpetraHelpers::BlockedTpetraOperator`. It takes an explicit list of monolithic global
IDs for each Teko block:

```cpp
#include "Teko_BlockedTpetraOperator.hpp"

std::vector<std::vector<Teko::GO>> blockGIDs(3);
const auto baseMap = A->getDomainMap();

for (size_t lid = 0; lid < baseMap->getLocalNumElements(); ++lid) {
  const Teko::GO gid = baseMap->getGlobalElement(static_cast<Teko::LO>(lid));
  if (gidBelongsToVelocity(gid)) {
    blockGIDs[0].push_back(gid);
  } else if (gidBelongsToPressure(gid)) {
    blockGIDs[1].push_back(gid);
  } else {
    blockGIDs[2].push_back(gid);
  }
}

auto blockedA = Teuchos::rcp(
    new Teko::TpetraHelpers::BlockedTpetraOperator(blockGIDs, A, "A"));
```

The outer vector entry is the Teko block number (`blockGIDs[0]` becomes block 0, and so on).
Each inner vector is this MPI rank's owned monolithic GID list for that block. The wrapped
operator is assumed to be square, with matching domain and range maps. The unit test
`tests/src/Tpetra/tBlockedTpetraOperator.cpp` is the current in-tree reference for constructing
these lists, checking `BlockedTpetraOperator::testAgainstFullOperator`, rebuilding after matrix
value changes, and applying reorder managers. A small runnable example built from that test
would be a useful addition to `examples/`.

---

## AddMultiplyPrecs — composing preconditioners

`examples/AddMultiplyPrecs/Driver.cpp` builds a block Gauss–Seidel factory and a block Jacobi
factory over the same diagonal strategy, then composes them:

```cpp
RCP<Teko::InverseFactory> inverse =
    Teko::invFactoryFromParamList(*paramList, "Amesos2");
RCP<Teko::BlockInvDiagonalStrategy> strategy =
    rcp(new Teko::InvFactoryDiagStrategy(inverse));

RCP<Teko::BlockPreconditionerFactory> GSFactory =
    rcp(new Teko::GaussSeidelPreconditionerFactory(Teko::GS_UseLowerTriangle, strategy));
RCP<Teko::BlockPreconditionerFactory> JacobiFactory =
    rcp(new Teko::JacobiPreconditionerFactory(strategy));

// Product of the two (MultPreconditionerFactory); AddPreconditionerFactory forms the sum
RCP<Teko::BlockPreconditionerFactory> MasterFactory =
    rcp(new Teko::MultPreconditionerFactory(GSFactory, JacobiFactory));
```

The equivalent parameter-driven form uses the `"Block Add"` / `"Block Multiply"`
[block types](04-block-preconditioners.md) with `"Preconditioner A"` / `"Preconditioner B"`.

**Full source:**
\include Driver.cpp

---

## step1 — writing your own BlockPreconditionerFactory

`examples/BuildPreconditioner/step1/` is the tutorial for extending Teko. To create a custom
block preconditioner you subclass `Teko::BlockPreconditionerFactory` and override one method,
`buildPreconditionerOperator`, using the block-operator utilities in `Teko_Utilities.hpp`:

```cpp
Teko::LinearOp ExamplePreconditionerFactory::buildPreconditionerOperator(
    Teko::BlockedLinearOp& blockOp, Teko::BlockPreconditionerState& state) const {
  // 1. Pull out the sub-blocks
  const Teko::LinearOp A_00 = Teko::getBlock(0, 0, blockOp);
  const Teko::LinearOp A_01 = Teko::getBlock(0, 1, blockOp);
  const Teko::LinearOp A_10 = Teko::getBlock(1, 0, blockOp);
  const Teko::LinearOp A_11 = Teko::getBlock(1, 1, blockOp);

  // 2. Cheap inverse of the (1,1) block, and an explicit (0,0) approximation
  const Teko::LinearOp invH = Teko::getInvDiagonalOp(A_11);
  const Teko::LinearOp P    = Teko::explicitAdd(A_00, Teko::scale(alpha_, A_01));
  const Teko::LinearOp invP = Teko::buildInverse(*inverse_, P);  // uses the injected inverse

  // 3. Assemble a block lower-triangular inverse operator
  Teko::BlockedLinearOp L = Teko::zeroBlockedOp(blockOp);
  Teko::setBlock(1, 0, L, A_10);
  Teko::endBlockFill(L);

  std::vector<Teko::LinearOp> invDiag = {invP, invH};
  return Teko::createBlockLowerTriInverseOp(L, invDiag);
}
```

Key utilities on display (all in `Teko_Utilities.hpp`): `blockRowCount`/`blockColCount`,
`getBlock`, `getInvDiagonalOp`, `explicitAdd`, `scale`, `buildInverse`, `zeroBlockedOp` /
`setBlock` / `endBlockFill`, and `createBlockLowerTriInverseOp`. The companion
`example-test.cpp` shows how to assemble four Tpetra sub-blocks, wrap them with
`Thyra::tpetraLinearOp`, combine with `Teko::block2x2`, and apply the resulting
preconditioner. See [Advanced Topics](07-advanced.md) for how to register such a factory
so it can be selected by a `"Type"` string, and for the state/rebuild pattern that extends
this example to cache inverse and explicit operators across nonlinear or time-step updates.

**Full source** (the custom factory):
\include ExamplePreconditionerFactory.cpp

**Full source** (the driver that exercises it):
\include example-test.cpp
