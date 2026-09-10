# Getting Started {#teko_getting_started}

[TOC]

This page introduces the core objects and walks through the minimal path to a working
block preconditioner: **assemble a blocked operator → choose a strategy → build the inverse
→ apply it**. It follows `examples/BuildPreconditioner/example-driver-belos.cpp` (the full,
runnable version).

## Core objects

Teko is layered on top of [Thyra](https://trilinos.github.io/thyra.html). Its central
operator and vector types are `Teuchos::RCP` handles around Thyra abstractions:

| Teko type | Underlying Thyra object held by the `RCP` | Role |
|-----------|------------------------------------------|------|
| `Teko::LinearOp` | `const Thyra::LinearOpBase<double>` | Any linear operator — a single block or a whole blocked system. |
| `Teko::BlockedLinearOp` | `Thyra::PhysicallyBlockedLinearOpBase<double>` | A block operator you assemble block-by-block. |
| `Teko::MultiVector` | `Thyra::MultiVectorBase<double>` | Solution/right-hand-side vectors, possibly blocked. |
| `Teko::BlockedMultiVector` | `Thyra::ProductMultiVectorBase<double>` | A product multi-vector with one component per block. |
| `Teko::VectorSpace` | `const Thyra::VectorSpaceBase<double>` | Domain/range space metadata for operators and vectors. |
| `Teko::InverseFactory` | — | Produces an operator approximating the inverse of a `LinearOp`. |
| `Teko::InverseLibrary` | — | A registry of *named* `InverseFactory` objects, built from a parameter list. |

Because Teko works on the Thyra layer, you wrap your native Tpetra matrices with
`Thyra::tpetraLinearOp` once, then work exclusively with `Teko::LinearOp`.

Two factory patterns are central to Teko. A `BlockPreconditionerFactory` knows how to build a
block preconditioner for a particular blocked operator, but it is reusable: the same factory can
be applied to different operators with the same mathematical structure. An `InverseFactory`
plays the same role for sub-solves. Block methods such as Jacobi, Gauss-Seidel, SIMPLE, and LSC
ask inverse factories to approximate the inverse of diagonal blocks or Schur-complement
operators without hard-coding whether that inverse is a direct solve, Krylov solve, single-level
preconditioner, or multigrid preconditioner.

## Step 1 — Build a blocked operator

Wrap your native Tpetra matrices as Thyra operators, then assemble them into a
block operator. The simplest constructor is `Thyra::block2x2`:

```cpp
// Mat and zeroMat are Tpetra::Operator wrapped as Thyra ops
Teko::LinearOp thMat  = Thyra::tpetraLinearOp<double>(range, domain, Mat);
Teko::LinearOp thZero = Thyra::tpetraLinearOp<double>(range, domain, zeroMat);

// A 2x2 block operator [ thMat  thZero ; thZero  thMat ]
Teko::LinearOp A = Thyra::block2x2(thMat, thZero, thZero, thMat);
```

For an arbitrary grid of blocks, use `Teko::zeroBlockedOp` / `Teko::setBlock` /
`Teko::endBlockFill` (see [Advanced Topics](07-advanced.md)). If instead you start from a
single *monolithic* matrix, let Teko split it for you: use **strided blocking** for regular
interleaved fields, or `BlockedTpetraOperator` for arbitrary field-to-GID maps (see below and
[Examples](06-examples.md)).

## Step 2 — Define named inverses

An `InverseLibrary` maps names to inverse factories. Build it from a parameter list:

```cpp
RCP<Teuchos::ParameterList> buildLibPL() {
  RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());

  Teuchos::ParameterList& sub_jac = pl->sublist("Jacobi");
  sub_jac.set("Type", "Block Jacobi");     // a Teko block preconditioner
  sub_jac.set("Inverse Type", "Ifpack2");  // each diagonal block solved with Ifpack2

  Teuchos::ParameterList& sub_gs = pl->sublist("Gauss-Seidel");
  sub_gs.set("Type", "Block Gauss-Seidel");
  sub_gs.set("Use Upper Triangle", true);
  sub_gs.set("Inverse Type", "Ifpack2");

  return pl;
}
```

Each sublist name (`"Jacobi"`, `"Gauss-Seidel"`) is a label *you* choose. The mandatory
`"Type"` key selects the backend. Here `"Ifpack2"` refers to a Stratimikos preconditioner;
because the library is built with a Stratimikos builder, all Stratimikos solvers and
preconditioners are available by name automatically. The full rules are on
[The Inverse Library](03-inverse-library.md) page.

```cpp
RCP<Stratimikos::DefaultLinearSolverBuilder> linearSolverBuilder =
    rcp(new Stratimikos::DefaultLinearSolverBuilder);

RCP<Teko::InverseLibrary> invLib =
    Teko::InverseLibrary::buildFromParameterList(*buildLibPL(), linearSolverBuilder);

RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("Gauss-Seidel");
```

## Step 3 — Build the inverse operator

`Teko::buildInverse` turns the factory + the blocked operator into a single `Teko::LinearOp`
that approximates \f$A^{-1}\f$:

```cpp
Teko::LinearOp prec = Teko::buildInverse(*inverse, A);
```

That `prec` is an ordinary Thyra operator. There is no separate "apply" ceremony — applying
`prec` *is* applying the preconditioner.

## Step 4 — Use it in a Krylov solve

Because `prec` and `A` are Thyra operators, they drop straight into Belos:

```cpp
RCP<Belos::LinearProblem<double, MV, OP> > problem =
    rcp(new Belos::LinearProblem<double, MV, OP>(A, x, b));
problem->setLeftPrec(prec);      // or setRightPrec(prec)
problem->setProblem();

RCP<Belos::SolverManager<double, MV, OP> > solver =
    rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcpFromRef(belosList)));
solver->solve();
```

## The native (non-Thyra) path

If you would rather keep your solver in Tpetra land, Teko provides operator wrappers
that build and apply the preconditioner without exposing Thyra:

```cpp
// vec describes the strided blocking: 2 velocity + 1 pressure unknown per node
std::vector<int> vec = {2, 1};
RCP<Teko::TpetraHelpers::StridedTpetraOperator> sA =
    rcp(new Teko::TpetraHelpers::StridedTpetraOperator(vec, A));  // A is a Tpetra::CrsMatrix

RCP<Teko::InverseLibrary> invLib   = Teko::InverseLibrary::buildFromStratimikos();
RCP<Teko::InverseFactory> inverse  = invLib->getInverseFactory("Amesos2");
RCP<Teko::NS::LSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(inverse, true));
RCP<Teko::BlockPreconditionerFactory> precFact =
    rcp(new Teko::NS::LSCPreconditionerFactory(strategy));

Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
prec.buildPreconditioner(sA);    // prec is now a Tpetra::Operator you can hand to Belos
```

This is the pattern in `examples/BuildPreconditioner/example-driver.cpp`.

### Native Tpetra with arbitrary block membership

Use `StridedTpetraOperator` when a fixed per-node stride describes the fields. If your
monolithic matrix uses another ordering — for example, application-owned GIDs for velocity,
pressure, temperature, or material fields are not regularly interleaved — use
`Teko::TpetraHelpers::BlockedTpetraOperator` instead. Its constructor takes a
`std::vector<std::vector<GO>>`: the outer vector is the Teko block number, and each inner
vector lists the monolithic global IDs owned by this MPI rank that belong to that block.

```cpp
#include "Teko_BlockedTpetraOperator.hpp"

using GO = Teko::GO;
using LO = Teko::LO;

std::vector<std::vector<GO>> blockGIDs(3);  // block 0, block 1, block 2
const auto baseMap = A->getDomainMap();     // for square systems, also the range map

for (size_t localRow = 0; localRow < baseMap->getLocalNumElements(); ++localRow) {
  const GO gid = baseMap->getGlobalElement(static_cast<LO>(localRow));

  if (isFluidVelocity(gid)) {
    blockGIDs[0].push_back(gid);  // Teko block 0: velocity rows/columns
  } else if (isFluidPressure(gid)) {
    blockGIDs[1].push_back(gid);  // Teko block 1: pressure rows/columns
  } else if (isTemperature(gid)) {
    blockGIDs[2].push_back(gid);  // Teko block 2: temperature rows/columns
  }
}

auto blockedA = Teuchos::rcp(
    new Teko::TpetraHelpers::BlockedTpetraOperator(blockGIDs, A, "A"));

Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
prec.buildPreconditioner(blockedA);
```

Across all MPI ranks, the inner vectors should partition the monolithic domain/range map for
ordinary non-overlapping field splits. `BlockedTpetraOperator` assumes a square operator with
matching domain and range maps. The GIDs in an inner vector do **not** need to be contiguous or
strided; Teko builds contiguous block-local maps internally for each block. The wrapper remains
a `Tpetra::Operator`, so it can be passed to `TpetraBlockPreconditioner::buildPreconditioner`
and then to Belos in the same way as the strided wrapper. Useful diagnostics include
`blockedA->GetBlock(i,j)`, `blockedA->WriteBlocks(prefix)`,
`blockedA->testAgainstFullOperator(count,tol)`, `blockedA->Reorder(manager)`, and
`blockedA->RebuildOps()` after matrix values change with the same block structure.

## Where to go next

- You configured the preconditioner in C++ above. To drive the *same* choices from an XML
  file (the usual production path), read the [Configuration Model](02-configuration-model.md).
- To understand the `"Inverse Type"` naming and how it reaches Ifpack2/MueLu/Amesos2/Belos,
  read [The Inverse Library](03-inverse-library.md).
- For every knob of every preconditioner, see the
  [Block Preconditioner Reference](04-block-preconditioners.md).
