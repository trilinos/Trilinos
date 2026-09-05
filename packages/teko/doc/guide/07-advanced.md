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

## Reusing inverses and explicit operators across iterations

In a Newton / nonlinear loop the operator changes every iteration but its sparsity and the
preconditioner *structure* often do not. Rebuild the inverse in place instead of reconstructing
it, which lets factories reuse symbolic factorizations, prolongation operators, and cached
explicit matrices when the chosen backend supports reuse:

```cpp
Teko::InverseLinearOp invA = Teko::buildInverse(*inverse, A0);   // first iteration
// ... new operator A1 with the same structure ...
Teko::rebuildInverse(*inverse, A1, invA);                        // refresh in place
```

### Store per-operator state in `BlockPreconditionerState`

A `BlockPreconditionerFactory` may build preconditioners for many different blocked operators.
Do not store operator-specific objects, such as cached inverses, directly in the factory unless
the factory is intentionally single-use. Use the `BlockPreconditionerState& state` argument to
`buildPreconditionerOperator` instead; Teko pairs that state object with the specific operator
instance being built or rebuilt.

For cached inverse operators, a common pattern is:

```cpp
Teko::LinearOp ExamplePreconditionerFactory::buildPreconditionerOperator(
    Teko::BlockedLinearOp& blockOp, Teko::BlockPreconditionerState& state) const {
  const Teko::LinearOp A00 = Teko::getBlock(0, 0, blockOp);
  const Teko::LinearOp A01 = Teko::getBlock(0, 1, blockOp);

  const Teko::LinearOp P = Teko::explicitAdd(A00, Teko::scale(alpha_, A01));

  Teko::ModifiableLinearOp& invP = state.getModifiableOp("invP");
  if (invP == Teuchos::null) {
    invP = Teko::buildInverse(*inverse_, P);
  } else {
    Teko::rebuildInverse(*inverse_, P, invP);
  }

  // Use invP when assembling the returned preconditioner operator.
  // ...
}
```

`state.getModifiableOp("invP")` returns a reference to a stored `RCP`. Assigning through that
reference updates the object kept by the state. The next rebuild for the same preconditioner
instance retrieves the same handle and can refresh it in place.

If the default state is not enough, derive from `Teko::BlockPreconditionerState` and override
your factory's `buildPreconditionerState()` method to return the richer state object. Inside
`buildPreconditionerOperator`, dynamically cast the incoming state to that derived type before
using the extra fields.

### Reuse explicit operators when possible

The explicit utilities have overloads that take a destination `ModifiableLinearOp`, for example
`explicitAdd(A, B, dest)` and `explicitMultiply(A, B, dest)`. Use these when a preconditioner
forms the same auxiliary operator every nonlinear iteration and only the values change:

```cpp
Teko::ModifiableLinearOp& P = state.getModifiableOp("P");
P = Teko::explicitAdd(A00, Teko::scale(alpha_, A01), P);
```

If `P` is null, Teko allocates it. If it is compatible with the new inputs, Teko can refill and
reuse it; otherwise it falls back to constructing a new operator. This is especially useful for
Schur-complement approximations and mass-scaled operators that would otherwise be reallocated
every rebuild.

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

## Testing and diagnostics

After implementing a new preconditioner factory, first test its action on known vectors before
using it inside a Krylov solve. A good unit test builds a small block operator, applies the
preconditioner to a predetermined vector, and compares the result with an independently computed
reference result (for example from a small script or a direct hand calculation). Then add an
integration test that uses the preconditioner inside Belos on a representative matrix.

Helpful diagnostics include:

- **`"Diagnostic Inverse"`** ([reference](04-block-preconditioners.md)) wraps any inner inverse
  to time it and optionally print its residual — drop it in temporarily to find the expensive
  sub-solve.
- **`"Test Block Operator"`** and **`"Write Block Operator"`**
  ([top-level parameters](02-configuration-model.md)) verify and dump the segregated blocks when
  you suspect the strided blocking or reordering is wrong.
- `StridedTpetraOperator::testAgainstFullOperator` and
  `BlockedTpetraOperator::testAgainstFullOperator` compare the wrapped block operator with the
  original monolithic `Tpetra::Operator` on random vectors.
- `Teuchos::describe(*op, Teuchos::VERB_EXTREME)` prints the structure and, for some concrete
  operators, values of a `Teko::LinearOp`.

## Thyra / Teko interoperability notes

Teko's core operator and vector types are aliases for `Teuchos::RCP` smart-pointer handles
around Thyra objects:

```cpp
Teko::VectorSpace        // Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
Teko::MultiVector        // Teuchos::RCP<Thyra::MultiVectorBase<double> >
Teko::BlockedMultiVector // Teuchos::RCP<Thyra::ProductMultiVectorBase<double> >
Teko::LinearOp           // Teuchos::RCP<const Thyra::LinearOpBase<double> >
Teko::BlockedLinearOp    // Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> >
```

A `BlockedLinearOp` handle can be upcast to a `LinearOp` handle, and a `BlockedMultiVector`
handle can be upcast to a `MultiVector` handle. The helper casts in `Teko_Utilities.hpp` make
those relationships explicit:

```cpp
Teko::LinearOp        A  = Teko::toLinearOp(blockedA);        // upcast; always valid
Teko::BlockedLinearOp BA = Teko::toBlockedLinearOp(A);        // downcast; throws if not blocked
Teko::MultiVector     x  = Teko::toMultiVector(blockedX);     // upcast; always valid
Teko::BlockedMultiVector BX = Teko::toBlockedMultiVector(x);  // downcast; throws if not blocked
```

Because these types are `Teuchos::RCP` handles, assignment is a shallow copy: it copies the
handle, not the matrix or vector data. Use the explicit `deepcopy` helpers when you need an
independent vector object:

```cpp
Teko::MultiVector copy = Teko::deepcopy(x);
Teko::BlockedMultiVector blockCopy = Teko::deepcopy(blockedX);
```

A few small utilities are often useful when constructing custom operators:

```cpp
Teko::LinearOp F = buildFOperator();
Teko::LinearOp I = Teko::identity(Teko::rangeSpace(F));

std::cout << Teuchos::describe(*F, Teuchos::VERB_EXTREME) << std::endl;
```

Wrap native matrices exactly once (`Thyra::tpetraLinearOp`), then stay in Teko/Thyra types.
Block assembly (`block2x2`, `zeroBlockedOp`/`setBlock`/`endBlockFill`) produces
`Teko::BlockedLinearOp`, which the factories consume directly. The `TpetraBlockPreconditioner`
wrapper exists so a finished Teko preconditioner can be handed back to a solver that only speaks
the Tpetra `Operator` interface.
