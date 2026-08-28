# Advanced Topics {#teko_advanced}

[TOC]

This page collects the mechanisms you reach for once the built-in preconditioners and their
options are not quite enough: supplying operators through the RequestHandler, reordering
blocks, reusing expensive inverses across nonlinear iterations, and registering your own
preconditioner factory so it can be named from XML.

## The RequestHandler / RequestCallback system

Some preconditioners need operators that cannot be derived from the assembled system — for
example a velocity mass matrix (SIMPLE/LSC), a pressure Laplacian, or a probing graph. Rather
than complicate every factory's constructor, Teko uses a **request/callback** pattern:

- A factory issues a `Teko::RequestMesg` (a named request, e.g. `"Velocity Mass Matrix"`)
  when it needs the operator.
- The application registers a `Teko::RequestCallback<DataT>` with a `RequestHandler` that
  answers matching messages.

The callback base (`src/Teko_RequestCallback.hpp`) has three methods:

```cpp
template <typename DataT>
class RequestCallback : public RequestCallbackBase {
 public:
  virtual bool  handlesRequest(const RequestMesg&) = 0;  // can I answer this message?
  virtual DataT request(const RequestMesg&)        = 0;  // produce the data
  virtual void  preRequest(const RequestMesg&)     = 0;  // prepare ahead of time (override; empty body if unused)
};
```

A minimal callback that supplies a velocity mass matrix:

```cpp
class MassMatrixCallback : public Teko::RequestCallback<Teko::LinearOp> {
 public:
  explicit MassMatrixCallback(Teko::LinearOp M) : massOp_(M) {}
  bool handlesRequest(const Teko::RequestMesg& m) override {
    return m.getName() == "Velocity Mass Matrix";
  }
  Teko::LinearOp request(const Teko::RequestMesg&) override { return massOp_; }
  void preRequest(const Teko::RequestMesg&) override {}
 private:
  Teko::LinearOp massOp_;
};
```

Register it and attach the handler to whatever builds the preconditioner:

```cpp
RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);
rh->addRequestCallback(rcp(new MassMatrixCallback(massOp)));

// Inside Stratimikos:
stratimikosFactory.setRequestHandler(rh);   // Teko::StratimikosFactory
// Or directly on a library:
invLib->setRequestHandler(rh);               // Teko::InverseLibrary
```

### Request message names used by the built-ins

| Message | Issued by | Enabled when |
|---------|-----------|--------------|
| `"Velocity Mass Matrix"` | LSC Basic Inverse, SIMPLE | `"Use Mass Scaling" = true` |
| `"W-Scaling Vector"` | LSC Basic Inverse | `"Use W-Scaling" = true` |
| `"Pressure Laplace Operator"` | LSC Pressure Laplace, PCD | always (for those strategies) |
| `"Velocity Mass Operator"` | LSC Pressure Laplace | always |
| `"PCD Operator"` | PCD strategy | always |
| `"Probing Graph"` | Probing Preconditioner | `"User Will Set Probing Graph" = true` |

The `"Required Parameters"` sublist of a Stratimikos-preconditioner
[library entry](03-inverse-library.md) is the parameter-driven counterpart: its contents
are removed and re-supplied through the RequestHandler at build time.

## Block reordering

Beyond the top-level `"Reorder Type"` string ([Configuration Model](02-configuration-model.md)),
you can reorder a block operator programmatically. `Teko::blockedReorderFromString` parses the
nested-bracket grammar into a `BlockReorderManager`, and `Teko::buildReorderedLinearOp`
applies it (`src/Teko_BlockedReordering.hpp`):

```cpp
std::string desc = "[2 [0 1]]";
RCP<const Teko::BlockReorderManager> mgr = Teko::blockedReorderFromString(desc);
Teko::LinearOp reordered = Teko::buildReorderedLinearOp(*mgr, blockedA);
```

The bracket string mirrors the block tree: bare integers are leaf block indices, and
`[ ... ]` groups its contents into a nested block. `"[2 [0 1]]"` builds a 2×2 outer operator
whose (0,0) entry is block 2 and whose (1,1) entry is the nested 2×2 of blocks 0 and 1.

## Reusing inverses across iterations

In a Newton / nonlinear loop the operator changes every iteration but its sparsity and the
preconditioner *structure* do not. Rebuild the inverse in place instead of reconstructing it,
which lets factories reuse symbolic factorizations and cached explicit operators:

```cpp
Teko::InverseLinearOp invA = Teko::buildInverse(*inverse, A0);   // first iteration
// ... new operator A1 with the same structure ...
Teko::rebuildInverse(*inverse, A1, invA);                        // refresh in place
```

`BlockPreconditionerState` (passed to `buildPreconditionerOperator`) is the per-instance
scratch space where a custom factory can stash explicit operators it wants to reuse between
rebuilds.

## Registering a custom preconditioner factory

Once you have written a `BlockPreconditionerFactory`
(see the [step1 example](06-examples.md)), register it so it can be selected by a
`"Type"` string from the [Inverse Library](03-inverse-library.md):

```cpp
Teko::PreconditionerFactory::addPreconditionerFactory(
    "My Block Prec",
    rcp(new Teko::AutoClone<MyBlockPreconditionerFactory>()));
```

After this, `"Type" = "My Block Prec"` works anywhere a Teko block type is expected, and your
factory's `initializeFromParameterList` receives the entry's remaining keys. Override
`initializeFromParameterList` to read your own parameters, and (optionally) provide
`getRequestHandler()`-based lookups for any operators you need from the application.

## Diagnostics

- **`"Diagnostic Inverse"`** ([reference](04-block-preconditioners.md)) wraps any inner inverse
  to time it and optionally print its residual — drop it in temporarily to find the expensive
  sub-solve.
- **`"Test Block Operator"`** and **`"Write Block Operator"`**
  ([top-level parameters](02-configuration-model.md)) verify and dump the segregated blocks when
  you suspect the strided blocking or reordering is wrong.

## Thyra / Teko interoperability notes

- `Teko::LinearOp` *is* a `const Thyra::LinearOpBase<double>` handle; any Thyra operator is a
  Teko operator and vice versa.
- Wrap native matrices exactly once (`Thyra::tpetraLinearOp`), then stay in Teko/Thyra types.
- Block assembly (`block2x2`, `zeroBlockedOp`/`setBlock`/`endBlockFill`) produces
  `Teko::BlockedLinearOp`, which the factories consume directly.
- The `TpetraBlockPreconditioner` wrapper exists so a finished Teko preconditioner can be
  handed back to a solver that only speaks the Tpetra `Operator` interface.
