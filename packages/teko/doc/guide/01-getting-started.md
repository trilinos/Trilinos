# Getting Started {#teko_getting_started}

[TOC]

This page introduces the core objects and walks through the minimal path to a working
block preconditioner: **assemble a blocked operator → choose a strategy → build the inverse
→ apply it**. It follows `examples/BuildPreconditioner/example-driver-belos.cpp` (the full,
runnable version).

## Core objects

Teko is layered on top of [Thyra](https://trilinos.github.io/thyra.html). Its central types
are typedefs into the Thyra abstractions:

| Teko type | Thyra type | Role |
|-----------|-----------|------|
| `Teko::LinearOp` | `Thyra::LinearOpBase<double>` (const) | Any linear operator — a single block or a whole blocked system. |
| `Teko::BlockedLinearOp` | `Thyra::PhysicallyBlockedLinearOpBase<double>` | A block operator you assemble block-by-block. |
| `Teko::MultiVector` | `Thyra::MultiVectorBase<double>` | Solution/right-hand-side vectors, possibly blocked. |
| `Teko::InverseFactory` | — | Produces an operator approximating the inverse of a `LinearOp`. |
| `Teko::InverseLibrary` | — | A registry of *named* `InverseFactory` objects, built from a parameter list. |

Because Teko works on the Thyra layer, you wrap your native Tpetra matrices with
`Thyra::tpetraLinearOp` once, then work exclusively with `Teko::LinearOp`.

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
single *monolithic* matrix whose unknowns are interleaved per node, you do not assemble
blocks by hand — you let Teko split it for you using **strided blocking** (see the
[Configuration Model](02-configuration-model.md) and the
[strided-operator examples](06-examples.md)).

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

## Where to go next

- You configured the preconditioner in C++ above. To drive the *same* choices from an XML
  file (the usual production path), read the [Configuration Model](02-configuration-model.md).
- To understand the `"Inverse Type"` naming and how it reaches Ifpack2/MueLu/Amesos2/Belos,
  read [The Inverse Library](03-inverse-library.md).
- For every knob of every preconditioner, see the
  [Block Preconditioner Reference](04-block-preconditioners.md).
