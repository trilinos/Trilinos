# The Inverse Library {#teko_inverse_library}

[TOC]

The **Inverse Library** (`Teko::InverseLibrary`) is the naming and dispatch engine that lets
one line of configuration — `"Inverse Type" = "MyThing"` — resolve to a direct solve, a
Krylov solver, an algebraic multigrid preconditioner, or another Teko block preconditioner.
Understanding it is what makes the rest of Teko's options composable.

## What a library entry is

Every sublist inside `"Inverse Factory Library"` is one entry. Its **name is a label you
choose**, and it must contain a `"Type"` key that names a backend
(`InverseLibrary::addInverse`, `src/Teko_InverseLibrary.cpp`). Everything else in the entry
is forwarded to that backend.

```xml
<ParameterList name="MyILU">          <!-- label: "MyILU" (your choice)          -->
  <Parameter name="Type" value="Ifpack2"/>   <!-- backend selector (required)      -->
  <Parameter name="Prec Type" value="ILUT"/> <!-- forwarded to Ifpack2             -->
</ParameterList>
```

You then reference the label wherever an inverse is expected: at the top level
(`"Inverse Type" = "MyILU"`) or inside a block preconditioner
(`"Inverse Type 1" = "MyILU"`).

## The three backend classes

`"Type"` is matched, in this order, against three internal lists (seeded in the
`InverseLibrary` constructor and extended with whatever the live Stratimikos builder
reports):

| Class | Seeded names | What Teko builds | Extra settings are… |
|-------|-------------|------------------|---------------------|
| **Stratimikos preconditioner** | `ML`, `Ifpack`, `Ifpack2`, `MueLu`, `Neumann Series` | a `PreconditionerInverseFactory` wrapping a `Thyra::PreconditionerFactoryBase` | that backend's own preconditioner parameters |
| **Stratimikos solver** | `Belos`, `Amesos`, `Amesos2`, `AztecOO` | a `SolveInverseFactory` (a full linear solve used as the inverse) | that solver's parameters (plus `"Use Preconditioner"`, below) |
| **Teko block preconditioner** | `Block Jacobi`, `Block Gauss-Seidel`, `NS SIMPLE`, `NS LSC`, `Block LU2x2`, … | the corresponding block factory | that factory's options ([reference](04-block-preconditioners.md)) |

If a `"Type"` matches none of the three, or a referenced label does not exist,
`getInverseFactory` prints the list of available names and aborts — a useful way to discover
what is registered.

## Bootstrapping: `buildFromStratimikos` vs `buildFromParameterList`

There are two ways to obtain a library:

- **`Teko::InverseLibrary::buildFromStratimikos()`** — auto-registers *every* Stratimikos
  solver and preconditioner under its own name. After this you can immediately ask for
  `getInverseFactory("Amesos2")`, `getInverseFactory("Ifpack2")`, etc. with no library
  entries at all. This is what the strided-operator examples use.

- **`Teko::InverseLibrary::buildFromParameterList(pl, builder)`** — builds the named entries
  from `pl` (the `"Inverse Factory Library"` sublist). Because you pass the Stratimikos
  builder, the built-in Stratimikos names remain available *in addition to* your entries, so
  entries can reference `"Ifpack2"` / `"Amesos2"` directly without defining them.

When Teko runs inside Stratimikos (the XML path), `StratimikosFactory` calls
`buildFromParameterList` on the `"Inverse Factory Library"` sublist for you — you never
construct the library by hand.

## How references chain

Any parameter whose *value* is a label is resolved by calling `getInverseFactory(label)` on
the same library. These parameters include, among others:

`"Inverse Type"`, `"Inverse Type <N>"`, `"Inverse Velocity Type"`, `"Inverse Pressure Type"`,
`"Inverse F Type"`, `"Inverse Laplace Type"`, `"Preconditioner Type"`,
`"Preconditioner Type <N>"`, `"Inverse Factory"`, `"Preconditioner A"`,
`"Preconditioner B"`, `"Use Preconditioner"`.

This is why a block preconditioner can point each of its sub-block solves at a different
entry — and those entries can themselves be block preconditioners, to any depth.

## Two special keys

### `"Use Preconditioner"` (solver entries)

Inside a Stratimikos **solver** entry (e.g. `Type = "Belos"`), the key `"Use Preconditioner"`
names *another* library entry to attach as that solver's preconditioner. This is how you
build a "preconditioned Krylov inverse" for a block:

```xml
<ParameterList name="Inverse Factory Library">
  <ParameterList name="BlockSolve">
    <Parameter name="Type" value="Belos"/>
    <Parameter name="Use Preconditioner" value="BlockPrec"/>   <!-- refers to entry below -->
    <ParameterList name="Solver Types"> <!-- Belos settings ... --> </ParameterList>
  </ParameterList>
  <ParameterList name="BlockPrec">
    <Parameter name="Type" value="Ifpack2"/>
    <Parameter name="Prec Type" value="RILUK"/>
  </ParameterList>
</ParameterList>
```

### `"Required Parameters"` (preconditioner entries)

Inside a Stratimikos **preconditioner** entry, a `"Required Parameters"` sublist is stripped
out before the backend sees it and re-supplied at build time through the
[RequestHandler](07-advanced.md). Use it to hand a backend operators/data that only the
application can provide (e.g. coordinates or a nullspace for MueLu).

## Worked example: nested solves (a Krylov-preconditioned block inverse)

Because any inverse-valued parameter can name another library entry, you can make each
sub-block of a block preconditioner be solved by a *full Krylov solver* that is itself
preconditioned — i.e. a genuine nested solve. The example below configures a
**Block Gauss-Seidel** whose two diagonal blocks are each solved by **GMRES (Belos)**, and
that inner GMRES is preconditioned by **Ifpack2 RILUK**:

```xml
<!-- Stratimikos "Preconditioner Types" -> "Teko" -->
<ParameterList name="Teko">
  <Parameter name="Inverse Type"     value="GS"/>
  <Parameter name="Strided Blocking" value="2 1"/>   <!-- e.g. 2 velocity + 1 pressure dof/node -->
  <ParameterList name="Inverse Factory Library">

    <!-- Outer block preconditioner: each diagonal block solved by a preconditioned GMRES -->
    <ParameterList name="GS">
      <Parameter name="Type"           value="Block Gauss-Seidel"/>
      <Parameter name="Inverse Type 1" value="PrecGMRES"/>   <!-- (0,0) block -->
      <Parameter name="Inverse Type 2" value="PrecGMRES"/>   <!-- (1,1) block -->
    </ParameterList>

    <!-- A full GMRES solve, used as a block inverse, preconditioned by "ILU" below -->
    <ParameterList name="PrecGMRES">
      <Parameter name="Type"               value="Belos"/>
      <Parameter name="Use Preconditioner" value="ILU"/>
      <Parameter name="Solver Type"        value="Block GMRES"/>
      <ParameterList name="Solver Types">
        <ParameterList name="Block GMRES">
          <Parameter name="Maximum Iterations"    value="50"/>
          <Parameter name="Convergence Tolerance" value="1.0e-4"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <!-- Preconditioner for the inner GMRES -->
    <ParameterList name="ILU">
      <Parameter name="Type"      value="Ifpack2"/>
      <Parameter name="Prec Type" value="RILUK"/>
    </ParameterList>

  </ParameterList>
</ParameterList>
```

The resolution chain is `GS` → (per block) `PrecGMRES` → (its preconditioner) `ILU`. To use a
*different* inner solver per block, point `"Inverse Type 1"` and `"Inverse Type 2"` at
distinct entries. The same pattern nests arbitrarily: a `PrecGMRES` entry could name another
block preconditioner as its `"Use Preconditioner"`, and so on.

## Building a library entry in C++

The XML/ParameterList path is preferred, but the same entries can be built directly. Two
convenience helpers exist:

```cpp
// From an already-populated "Inverse Factory Library" style parameter list:
RCP<Teko::InverseFactory> inv = Teko::invFactoryFromParamList(paramList, "Amesos2");

// Or grab a Stratimikos default by name:
RCP<Teko::InverseLibrary> lib  = Teko::InverseLibrary::buildFromStratimikos();
RCP<Teko::InverseFactory> inv2 = lib->getInverseFactory("Ifpack2");
```

Once you have an `InverseFactory`, `Teko::buildInverse(*inv, A)` produces the inverse
operator (and `Teko::rebuildInverse(*inv, A, invOp)` refreshes an existing one in place — see
[Advanced Topics](07-advanced.md) for reuse).
