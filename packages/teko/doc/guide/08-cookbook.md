# Cookbook {#teko_cookbook}

[TOC]

This page collects common Teko usage patterns. Use it as a starting point when deciding how
to block an operator, which inverse-library entries to define, and where an application must
supply extra physics operators.

## Choose the blocking route

| Starting point | Use | Notes |
|----------------|-----|-------|
| Separate submatrices already exist | `Thyra::block2x2` or `Teko::zeroBlockedOp` / `Teko::setBlock` / `Teko::endBlockFill` | Best when the application naturally assembles field-to-field blocks. |
| One monolithic Tpetra matrix with a fixed interleaved nodal layout | `Teko::TpetraHelpers::StridedTpetraOperator` or XML `"Strided Blocking"` | Example: `[u v p]` at each node gives `{2, 1}` in C++ or `"2 1"` in XML. |
| One monolithic Tpetra matrix with arbitrary field ownership/order | `Teko::TpetraHelpers::BlockedTpetraOperator` | Pass a `std::vector<std::vector<GO>>`; outer index is the Teko block number, inner entries are monolithic global IDs owned on this rank. |

### Arbitrary Tpetra blocks from global IDs

For non-strided application orderings, build a block-membership list from the monolithic map:

```cpp
#include "Teko_BlockedTpetraOperator.hpp"

std::vector<std::vector<Teko::GO>> blockGIDs(3);  // velocity, pressure, temperature
const auto baseMap = A->getDomainMap();          // for square systems, also the range map

for (size_t lid = 0; lid < baseMap->getLocalNumElements(); ++lid) {
  const Teko::GO gid = baseMap->getGlobalElement(static_cast<Teko::LO>(lid));

  if (isVelocityGID(gid)) {
    blockGIDs[0].push_back(gid);  // Teko block 0
  } else if (isPressureGID(gid)) {
    blockGIDs[1].push_back(gid);  // Teko block 1
  } else if (isTemperatureGID(gid)) {
    blockGIDs[2].push_back(gid);  // Teko block 2
  }
}

auto blockedA = Teuchos::rcp(
    new Teko::TpetraHelpers::BlockedTpetraOperator(blockGIDs, A, "A"));
```

Each rank lists only the GIDs it owns. The inner vectors should form the desired field split of
the monolithic map; for the usual non-overlapping case, a GID appears in exactly one inner
vector. `BlockedTpetraOperator` assumes a square operator with matching domain and range maps.
Teko constructs contiguous block-local maps internally, extracts every `(i,j)` subblock, and
keeps a Tpetra `Operator` interface for the original monolithic map.

Use it exactly like the strided wrapper when building a native Tpetra preconditioner:

```cpp
Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFactory);
prec.buildPreconditioner(blockedA);

RCP<Tpetra::Operator<Teko::ST, Teko::LO, Teko::GO, Teko::NT>> precOp =
    Teuchos::rcp(&prec, false);
belosProblem->setRightPrec(precOp);
```

Useful checks while developing the GID lists:

```cpp
blockedA->testAgainstFullOperator(20, 1.0e-12);
blockedA->WriteBlocks("A-block");
auto A01 = blockedA->GetBlock(0, 1);
```

## Generic two-field block systems

For a generic two-field block system with nonsingular diagonal blocks

\f[
\begin{bmatrix} A_{00} & A_{01} \\ A_{10} & A_{11} \end{bmatrix},
\f]

block Gauss-Seidel is often a useful starting point: define one inverse per diagonal block and
apply them in a triangular sweep.

```xml
<ParameterList name="Teko">
  <Parameter name="Strided Blocking" type="string" value="2 1"/>
  <Parameter name="Inverse Type"     type="string" value="BGS"/>
  <ParameterList name="Inverse Factory Library">
    <ParameterList name="BGS">
      <Parameter name="Type"           type="string" value="Block Gauss-Seidel"/>
      <Parameter name="Inverse Type 1" type="string" value="Block0Solve"/>
      <Parameter name="Inverse Type 2" type="string" value="Block1ILU"/>
    </ParameterList>
    <ParameterList name="Block0Solve">
      <Parameter name="Type" type="string" value="Belos"/>
      <Parameter name="Use Preconditioner" type="string" value="Block0MG"/>
    </ParameterList>
    <ParameterList name="Block0MG">
      <Parameter name="Type" type="string" value="MueLu"/>
    </ParameterList>
    <ParameterList name="Block1ILU">
      <Parameter name="Type" type="string" value="Ifpack2"/>
      <Parameter name="Prec Type" type="string" value="RILUK"/>
    </ParameterList>
  </ParameterList>
</ParameterList>
```

Use `Block Jacobi` if you need a cheaper or more parallel preconditioner.

## Saddle-point systems

For saddle-point systems, one diagonal block is often zero or singular, so the generic block
Jacobi/Gauss-Seidel recipe above is not the right model: it would require inverting that
singular diagonal block. Instead, use a Schur-complement or physics-specific approximation.
For a generic 2-by-2 saddle-point problem, start with `Block LU2x2` and provide an inverse for
the nonsingular field block plus an inverse for the Schur-complement approximation. For
incompressible Navier-Stokes systems, prefer the dedicated SIMPLE, LSC, or PCD strategies below.

## Navier-Stokes: SIMPLE, LSC, and PCD

For incompressible Navier-Stokes, prefer the built-in physics factories over manually writing a
saddle-point approximation:

- Use `"NS SIMPLE"` as a robust first configuration. It needs no extra operators unless
  `"Use Mass Scaling" = true`.
- Use `"NS LSC"` when a least-squares commutator approximation is appropriate. Some strategies
  request a velocity mass matrix, pressure Laplace operator, or W-scaling vector.
- Use `"Block LU2x2"` with `"Strategy Name" = "NS PCD Strategy"` when the application can
  provide the pressure operators required by PCD.

See [Navier-Stokes Preconditioners](05-navier-stokes.md) for the complete option tables and
request names.

## Flow plus energy, species, or other transported scalars

A common multiphysics pattern is a Navier-Stokes block coupled to one or more scalar transport
blocks. For a strided layout `[u v p T]`, first split into velocity, pressure, and temperature,
then use a reorder to group velocity-pressure into a nested Navier-Stokes subsystem:

```xml
<ParameterList name="Teko">
  <Parameter name="Strided Blocking" type="string" value="2 1 1"/>
  <Parameter name="Reorder Type"     type="string" value="[[0 1] 2]"/>
  <Parameter name="Inverse Type"     type="string" value="OuterGS"/>
  <ParameterList name="Inverse Factory Library">
    <ParameterList name="OuterGS">
      <Parameter name="Type"           type="string" value="Block Gauss-Seidel"/>
      <Parameter name="Inverse Type 1" type="string" value="FlowPrec"/>
      <Parameter name="Inverse Type 2" type="string" value="EnergyPrec"/>
    </ParameterList>
    <ParameterList name="FlowPrec">
      <Parameter name="Type" type="string" value="NS SIMPLE"/>
      <Parameter name="Inverse Velocity Type" type="string" value="VelocityMG"/>
      <Parameter name="Inverse Pressure Type" type="string" value="PressureILU"/>
    </ParameterList>
    <ParameterList name="EnergyPrec">
      <Parameter name="Type" type="string" value="Ifpack2"/>
      <Parameter name="Prec Type" type="string" value="RILUK"/>
    </ParameterList>
  </ParameterList>
</ParameterList>
```

After `"Strided Blocking" = "2 1 1"`, the flat block order is `0 = velocity`, `1 = pressure`,
`2 = temperature`. The reorder `"[[0 1] 2]"` makes a 2-by-2 outer operator whose block 0 is the
nested velocity-pressure operator and whose block 1 is the temperature operator. Therefore
`"Inverse Type 1"` in `OuterGS` applies to the nested flow subsystem, and `"Inverse Type 2"`
applies to the temperature block.

The same idea extends to arbitrary GID blocking: make `blockGIDs[0]` velocity, `blockGIDs[1]`
pressure, `blockGIDs[2]` temperature, construct a `BlockedTpetraOperator`, and call
`Reorder(*Teko::blockedReorderFromString("[[0 1] 2]"))` before building the preconditioner.

## Application-supplied operators

When a preconditioner needs physics information that is not present in the monolithic matrix,
register a `RequestHandler` callback. Typical requests are:

| Application data | Request message | Used by |
|------------------|-----------------|---------|
| Velocity mass matrix | `"Velocity Mass Matrix"` | SIMPLE and LSC mass scaling |
| W-scaling vector | `"W-Scaling Vector"` | LSC W-scaling |
| Pressure Laplace operator | `"Pressure Laplace Operator"` | LSC pressure-Laplace strategy and PCD |
| Velocity mass operator | `"Velocity Mass Operator"` | LSC pressure-Laplace strategy |
| PCD pressure convection-diffusion operator | `"PCD Operator"` | PCD strategy |
| Probing graph | `"Probing Graph"` | Probing preconditioner |

Keep callbacks small: they should return already-assembled operators or data owned by the
application, not rebuild expensive objects unless the nonlinear/time state requires it.

## Reuse in nonlinear or time-stepping loops

For repeated solves with the same block structure:

1. Keep the same factory and state object.
2. Update values in the application matrix.
3. If using a Tpetra blocking wrapper, call `RebuildOps()` so the extracted subblocks are
   refreshed.
4. Call `Teko::rebuildInverse(*factory, newA, oldInverse)` or
   `TpetraBlockPreconditioner::rebuildPreconditioner(...)`.
5. Cache any auxiliary explicit operators in `BlockPreconditionerState` using named
   `ModifiableLinearOp` entries.

This avoids repeated allocation and lets inner solver packages reuse symbolic setup when they
can.

## Debugging recipes

- If the preconditioner behaves like the wrong operator, enable `"Test Block Operator"` or call
  `testAgainstFullOperator` on a Tpetra blocking wrapper.
- If a field split is suspect, enable `"Write Block Operator"` in XML or call `WriteBlocks` and
  inspect the Matrix-Market files.
- If nested solves are unexpectedly expensive, wrap suspect entries with `"Diagnostic Inverse"`.
- If an XML entry is not found, intentionally check the error output from `getInverseFactory`;
  it prints the registered names from the active build.
- If a custom factory receives a non-blocked operator, use `Teko::toBlockedLinearOp` during
  development to fail early with a clear cast error.

## Documentation roadmap

A useful next step is to keep this guide organized around workflows rather than implementation
classes alone:

1. **Cookbook first.** Add one short recipe per application pattern, each with the blocking
   choice, a minimal XML or C++ snippet, the expected `RequestHandler` requests, and the
   diagnostics to run first.
2. **Reference second.** Keep exhaustive parameter tables in
   [Block Preconditioner Reference](04-block-preconditioners.md) and link to them from recipes
   instead of repeating every option.
3. **Examples as tests.** Whenever a recipe stabilizes, promote it to a small example under
   `packages/teko/examples/` so Doxygen can include real source and CI can catch API drift.
4. **Backend-neutral language.** Describe inner solvers as inverse-library labels and mention
   only currently supported Trilinos packages in examples.
5. **Reuse guidance near nonlinear examples.** Every transient or nonlinear recipe should state
   whether the matrix graph, field split, and auxiliary operators are reusable and which rebuild
   API to call.

## Candidate future examples

The documentation would benefit from these small, runnable additions under `packages/teko/examples/`: 

1. **ArbitraryBlockedTpetra** — construct `std::vector<std::vector<GO>>` from an application
   field map, wrap a monolithic matrix with `BlockedTpetraOperator`, verify with
   `testAgainstFullOperator`, and solve with Belos plus `TpetraBlockPreconditioner`.
2. **XMLBlockGaussSeidel** — a minimal Stratimikos XML file selecting Teko, splitting a
   two-field matrix, and assigning different inverse-library entries to the two diagonal
   blocks.
3. **RequestHandlerOperators** — a compact example that supplies a velocity mass matrix and a
   pressure Laplace operator to an LSC or PCD setup.
4. **ReuseLoop** — update matrix values in a pseudo time-step loop and demonstrate
   `RebuildOps()` plus `rebuildPreconditioner` / `rebuildInverse`.
5. **FlowEnergyNested** — demonstrate the flow-plus-temperature pattern above with a nested
   reorder and a physics preconditioner on the inner velocity-pressure block.

Keeping these examples small and matrix-file driven would make them useful both as Doxygen
snippets and as regression tests for the documented workflows.
