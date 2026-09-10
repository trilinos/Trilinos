# Configuration Model {#teko_config_model}

[TOC]

This is the most important page in the guide. Teko preconditioners are configured entirely
through nested `Teuchos::ParameterList`s (equivalently, Trilinos XML). Understanding the
three layers, and the one required top-level structure, is enough to configure any Teko
preconditioner.

## The three layers

```
Preconditioner Types
└── Teko                              (1) Stratimikos layer  — how to block the operator
    ├── "Strided Blocking" = "..."
    ├── "Reorder Type"     = "..."
    ├── "Inverse Type"     = "<name>"
    └── Inverse Factory Library       (2) the named-inverse registry
        ├── <name-A>
        │   ├── "Type" = "<backend>"
        │   └── ...                   (3) per-factory options
        └── <name-B>
            ├── "Type" = "<backend>"
            └── ...
```

1. **Stratimikos layer** — the `Preconditioner Types → Teko` sublist. Its keys describe how
   to break the monolithic operator into blocks and which named inverse to apply to the
   whole system.
2. **Inverse Factory Library** — a sublist that defines every *named* inverse. Each entry
   has a mandatory `"Type"` selecting a backend.
3. **Per-factory options** — the keys each Teko block preconditioner reads for itself.

## Layer 1 — top-level `Teko` parameters

These are the valid keys of the `Preconditioner Types → Teko` sublist. They are defined and
documented in `StratimikosFactory::getValidParameters()`
(`src/Teko_StratimikosFactory.cpp`), which means Stratimikos will validate them and reject
typos.

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `"Inverse Type"` | string | `"Amesos2"` if enabled, else `""` | Name of the inverse applied to the full blocked system. Either an available Stratimikos default name or a label defined in `"Inverse Factory Library"`. |
| `"Strided Blocking"` | string | `"1"` | How to split the operator into blocks. `"3 1 1"` groups the first 3 unknowns/node together, then two singletons — see below. |
| `"Reorder Type"` | string | `""` (no reorder) | Nested reordering of the blocks, e.g. `"[2 [0 1]]"`. See below. |
| `"Test Block Operator"` | bool | `false` | Diagnostic: apply the segregated operator and compare against the original. |
| `"Write Block Operator"` | bool | `false` | Diagnostic: write the segregated blocks to disk as `block-?_xx`. |
| `"Inverse Factory Library"` | sublist | — | Definitions of all named Teko inverses (layer 2). |

### Strided Blocking

`"Strided Blocking"` describes how the interleaved unknowns of a monolithic operator map to
blocks. The string is a space-separated list of block widths, "associated with the solution
vector." For 3D Navier–Stokes with an energy equation the per-node unknowns are
`[u v w p T]`; the string `"3 1 1"` means:

- put the first **3** unknowns per node (the velocity components `u v w`) into one block,
- put the next **1** (`p`) into its own block,
- put the next **1** (`T`) into its own block.

The result is a 3×3 block operator. `"1"` (the default) means every unknown is its own
block in a scalar sense — i.e. no strided grouping.

### Reorder Type

`"Reorder Type"` optionally nests the flat blocks into a hierarchy *after* strided blocking.
Given the `"3 1 1"` blocking above (blocks indexed `0`=velocity, `1`=pressure, `2`=temp),
the string `"[2 [0 1]]"` produces a 2×2 outer system whose (0,0) entry is itself the 2×2
velocity–pressure system and whose (1,1) entry is the temperature block. This lets you build
"a preconditioner of a preconditioner" — e.g. a Navier–Stokes method on the inner
velocity–pressure block wrapped in an outer block method that couples in temperature. The
grammar is the nested-bracket form parsed by `Teko::blockedReorderFromString`
(see [Advanced Topics](07-advanced.md)).

## Layer 2 — the Inverse Factory Library

Every entry of `"Inverse Factory Library"` is a sublist whose *name is a label you choose*
and whose mandatory `"Type"` selects a backend. Three kinds of backend exist:

- **Stratimikos solvers** — commonly `"Belos"` and `"Amesos2"` (a full linear solve as the
  "inverse"). The exact set is whatever the configured `Stratimikos::DefaultLinearSolverBuilder`
  reports as enabled.
- **Stratimikos preconditioners** — commonly `"Ifpack2"`, `"MueLu"`, and
  `"Neumann Series"` (a single-level/multilevel preconditioner as the "inverse"). The exact set
  is likewise taken from the active Stratimikos builder.
- **Teko block preconditioners** — `"Block Jacobi"`, `"Block Gauss-Seidel"`, `"NS SIMPLE"`,
  `"NS LSC"`, `"Block LU2x2"`, … (the block methods in the
  [reference](04-block-preconditioners.md)).

All remaining keys in an entry are passed to that backend. For a Stratimikos backend they
are the backend's own settings (e.g. Ifpack2's `"Prec Type"`); for a Teko block backend they
are that factory's options. The complete dispatch rules — and how one entry references
another — are on [The Inverse Library](03-inverse-library.md) page.

## A complete, minimal example

This is a full block Gauss–Seidel preconditioner driven entirely from parameters. It splits
a two-field operator (`"1 1"`), applies the named inverse `"BGS"`, and `"BGS"` in turn
solves each diagonal block with Amesos2. (Reconstructed from
`tests/unit_tests/tStratimikosFactory.cpp`, which builds exactly this list in C++.)

**As XML:**

```xml
<ParameterList name="Stratimikos">
  <Parameter name="Preconditioner Type" type="string" value="Teko"/>
  <ParameterList name="Preconditioner Types">
    <ParameterList name="Teko">
      <Parameter name="Write Block Operator" type="bool"   value="false"/>
      <Parameter name="Test Block Operator"  type="bool"   value="false"/>
      <Parameter name="Strided Blocking"     type="string" value="1 1"/>
      <Parameter name="Inverse Type"         type="string" value="BGS"/>
      <ParameterList name="Inverse Factory Library">
        <ParameterList name="BGS">
          <Parameter name="Type"         type="string" value="Block Gauss-Seidel"/>
          <Parameter name="Inverse Type" type="string" value="Amesos2"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>
</ParameterList>
```

**The equivalent C++ ParameterList:**

```cpp
Teuchos::ParameterList& tekoList =
    params->sublist("Preconditioner Types").sublist("Teko");
tekoList.set("Write Block Operator", false);
tekoList.set("Test Block Operator", false);
tekoList.set("Strided Blocking", "1 1");
tekoList.set("Inverse Type", "BGS");

Teuchos::ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
ifl.sublist("BGS").set("Type", "Block Gauss-Seidel");
ifl.sublist("BGS").set("Inverse Type", "Amesos2");
```

> **Note.** No XML file in the Teko source tree currently wires `"Teko"` as a
> `Preconditioner Type`; the two `examples/*/solverparams.xml` files configure only the
> inner Stratimikos solver. The example above is therefore the reference pattern to copy.

## Enabling Teko inside Stratimikos

For an application that already uses Stratimikos, register Teko once and then the XML above
just works:

```cpp
Stratimikos::DefaultLinearSolverBuilder stratBuilder;
Teko::addTekoToStratimikosBuilder(stratBuilder, "Teko");   // register the "Teko" preconditioner

stratBuilder.setParameterList(myParamsFromXml);            // the list shown above
RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
    stratBuilder.createPreconditioningStrategy("Teko");
```

`Teko::addTekoToStratimikosBuilder` is declared in `src/Teko_StratimikosFactory.hpp` and is
the single public entry point for making Teko available as a Stratimikos preconditioner.

## A richer example: nesting a block method over a multigrid block solve

The library can hold as many named inverses as you like, and block methods can reference
them. Here a block Gauss–Seidel preconditioner solves its (0,0) velocity block with MueLu
and its (1,1) pressure block with Ifpack2:

```xml
<ParameterList name="Teko">
  <Parameter name="Strided Blocking" type="string" value="2 1"/>
  <Parameter name="Inverse Type"     type="string" value="BGS"/>
  <ParameterList name="Inverse Factory Library">

    <ParameterList name="BGS">
      <Parameter name="Type"           type="string" value="Block Gauss-Seidel"/>
      <Parameter name="Use Upper Triangle" type="bool" value="false"/>
      <Parameter name="Inverse Type 1" type="string" value="VelocityMG"/>
      <Parameter name="Inverse Type 2" type="string" value="PressureILU"/>
    </ParameterList>

    <ParameterList name="VelocityMG">
      <Parameter name="Type" type="string" value="MueLu"/>
      <!-- MueLu's own settings go here -->
      <Parameter name="number of equations" type="int" value="2"/>
    </ParameterList>

    <ParameterList name="PressureILU">
      <Parameter name="Type"      type="string" value="Ifpack2"/>
      <Parameter name="Prec Type" type="string" value="ILUT"/>
      <ParameterList name="Ifpack2 Settings">
        <Parameter name="fact: ilut level-of-fill" type="double" value="1.0"/>
      </ParameterList>
    </ParameterList>

  </ParameterList>
</ParameterList>
```

`"Inverse Type 1"` / `"Inverse Type 2"` are the per-block overrides for block Gauss–Seidel
(the `<N>` is a 1-based block index). Their values (`"VelocityMG"`, `"PressureILU"`) are
themselves library labels, resolved recursively. This chaining is the heart of Teko's
flexibility and is documented next.
