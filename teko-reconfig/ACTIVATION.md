# How the Krylov-surrogate hack activates (and what breaks it)

No application API change. After a flexible, blocked GMRES solve,
`Belos::BlockGmresSolMgr::solve()` calls a hook that `libteko` registers at
load time (`Teko_KrylovSurrogateInit.cpp` → `Belos::AdaptiveHook::registerHook`).
The hook (`Teko::KrylovSurrogate::adaptiveLoop`) builds the surrogate, asks the
watcher for an ordering, and re-solves. A drop-in app replaces `teko_ext`
without calling anything new — it just has to do a solve that meets the
conditions below.

## Activates only when ALL hold
- **`TEKO_ADAPTIVE_RECONFIG` is set** (truthy). The hook is opt-in: even when
  registered it returns immediately unless this env var is set, so linking
  Teko does not silently hijack solves. The Python interface sets it.
- The app **links/loads `libteko`** (its static initializer registers the
  hook; otherwise `invoke()` is a silent no-op).
- Solve uses **`Belos::BlockGmresSolMgr`** with **`"Flexible Gmres": true`**.
- Template types are exactly **`SC=double`, `MV=Thyra::MultiVectorBase<double>`,
  `OP=Thyra::LinearOpBase<double>`** (the `if constexpr` guard in the sol mgr).
- The operator casts to **`Thyra::BlockedLinearOpBase<double>`** (physically
  blocked), and its blocks are **`Thyra::TpetraLinearOp`** over
  **`Tpetra::CrsMatrix<double,int,long long,OpenMP-node>`** (the aliases in
  `Teko_KrylovSurrogate.hpp`).
- The first solve produced Krylov data (`curDim > 0`).
- A **watcher** answers `s<N>_reconfig.json` in the requests dir
  (`TEKO_RECONFIG_REQUESTS_DIR`, else the hardcoded `kDefaultRequestsDir`).

## Breaks / silently does nothing when
- **`TEKO_ADAPTIVE_RECONFIG` not set** → hook returns immediately (by design).
- **`libteko` not loaded** → hook never registered → no-op. (Also: the
  registration is a file-scope static; a static-archive Teko would need
  `--whole-archive` or it gets stripped. Shared `libteko.so` is fine.)
- Wrong manager (`PseudoBlockGmres`, …) or `Flexible Gmres` false → no invoke.
- Non-`double`, or Tpetra-templated `MV`/`OP` → excluded by `if constexpr`.
- Operator not blocked, or blocks aren't `TpetraLinearOp`/`CrsMatrix` with
  matching `GO=long long`, `Node=OpenMP` → casts return null → those blocks are
  treated as zero → degenerate/wrong surrogate (no error raised).
- No watcher, or watcher writing to a different dir → 600 s timeout, then the
  first-solve result is returned unchanged. (The watcher is now managed by the
  interface: started once on rank 0, stopped at interpreter exit; as a safety
  net it self-terminates after 1000 s idle or if its parent dies.)
- `reconfig.json` lacking `"use_ordering"`, or a non-integer ordering entry →
  throws out of the solve.
- Running many solves per process leaks factor-time across solves unless the
  registry is reset at each solve start (`teko_ext` does this; a drop-in app
  must too — see `Teko::FactorTimeRegistry::reset()`).

Note: a candidate ordering whose factorization fails on only some ranks no
longer deadlocks the sweep — `solveOrdering` max-reduces a build-failure flag
and all ranks skip that candidate together (mirroring `teko_ext`'s
diagonal-block build).

## Optional knobs (env)
- `TEKO_ADAPTIVE_RECONFIG` — master opt-in; the hook is inert unless set.
- `TEKO_RECONFIG_REQUESTS_DIR` — handshake directory (must match the watcher).
- `TEKO_FACTOR_WARMUP` — if set, the first factorization runs once untimed so
  cold-start cost stays out of solve-1's reported factor time.
- `TEKO_WATCHER_IDLE_TIMEOUT` — watcher self-terminate-after-idle seconds
  (default 1000; ≤ 0 disables).
