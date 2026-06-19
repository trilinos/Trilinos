# How the Krylov-surrogate hack activates (and what breaks it)

The whole mechanism is gated on **`TEKO_ADAPTIVE_RECONFIG`** — set it to a
truthy value or nothing below happens (the hook returns immediately, so merely
linking Teko never hijacks a solve). With it set, after a flexible, blocked
GMRES solve `Belos::BlockGmresSolMgr::solve()` calls a hook that `libteko`
registers at load time (`Teko_KrylovSurrogateInit.cpp` →
`Belos::AdaptiveHook::registerHook`). The hook
(`Teko::KrylovSurrogate::adaptiveLoop`) builds the surrogate, gets a response
from the watcher (`use_ordering` plus optional `test_orderings` /
`selection_mode`), optionally builds-and-times each candidate ordering (written
to `s<N>_solved.json`), and re-solves with the selected one.

**This fires on both converged and stalled first solves.** A solve that hit its
max-iteration limit (`curDim > 0`) is exactly the case reconfiguration should
rescue, so it takes the same path: the surrogate builds fine from its Krylov
data, and the re-solve runs. If that re-solve converges, the hook reports it
back and `BlockGmresSolMgr::solve()` returns **converged** with the re-solve's
iteration count / residual — so a deliberately short, low-max-iter first solve
that stalls is rescued and the caller sees success rather than the original
stall. `s<N>_conv.json` records both `solve1` (the stall: iterations == max
iters) and `solve2` (the rescue). Caveat: the re-solve inherits the first
solve's parameters, including `Maximum Iterations` — so if max-iters is set
very low for cheap data-gathering, the re-solve is capped at that same low
value and may not converge; give it a larger budget if you need the rescue to
land. (A `curDim == 0` breakdown still returns early — nothing to build from.)

## Where the knobs are defined

| Env var | Read in | Default |
| --- | --- | --- |
| `TEKO_ADAPTIVE_RECONFIG` | [`Teko_KrylovSurrogate.hpp:1401`](../packages/teko/src/Teko_KrylovSurrogate.hpp#L1401) | unset → inert |
| `TEKO_RECONFIG_REQUESTS_DIR` | [`Teko_KrylovSurrogate.hpp:1409`](../packages/teko/src/Teko_KrylovSurrogate.hpp#L1409) | `kDefaultRequestsDir` ([`:1147`](../packages/teko/src/Teko_KrylovSurrogate.hpp#L1147)) |
| `TEKO_FACTOR_WARMUP` | [`Teko_InverseFactory.cpp:130`](../packages/teko/src/Teko_InverseFactory.cpp#L130) | on |
| `TEKO_WATCHER_IDLE_TIMEOUT` | [`wait_for_request.py:40`](wait_for_request.py#L40) | 1000 s |

## Activates only when ALL hold
- **`TEKO_ADAPTIVE_RECONFIG` is set** (truthy) — the gate above; inert otherwise.
- The app **links/loads `libteko`** (its static initializer registers the
  hook; otherwise `invoke()` is a silent no-op).
- Solve uses **`Belos::BlockGmresSolMgr`** with **`"Flexible Gmres": true`**.
- Template types are exactly **`SC=double`, `MV=Thyra::MultiVectorBase<double>`,
  `OP=Thyra::LinearOpBase<double>`** (the `if constexpr` guard in the sol mgr).
- The operator casts to **`Thyra::BlockedLinearOpBase<double>`** (physically
  blocked), and its blocks are **`Thyra::TpetraLinearOp`** over
  **`Tpetra::CrsMatrix<double,int,long long,OpenMP-node>`** (the aliases in
  `Teko_KrylovSurrogate.hpp`).
- The first solve produced Krylov data (`curDim > 0`) — true whether it
  converged or stalled at max-iters. A `curDim == 0` breakdown returns early.
- A **watcher** answers `s<N>_reconfig.json` in the requests dir
  (`TEKO_RECONFIG_REQUESTS_DIR`, else the hardcoded `kDefaultRequestsDir`).
  Needed for both the converged and the stalled (rescue) paths, since both
  build the surrogate and reconfigure.

## File I/O formats

All files live in the requests dir and are named `s<N>_*.json`, where `N` is a
monotonic id chosen by `nextRequestNumber()` = max existing id + 1 over
`s<N>_request.json` / `_conv.json` / `_solved.json` (so ids are never reused
within a run; the Python interface wipes them at startup so each run restarts
at `s0`). Every file is written to `*.tmp` and atomically renamed, so a reader
never sees a half-written file. Per request `N` the exchange is:

```
 C++  --s<N>_request.json-->  watcher
 C++  <--s<N>_reconfig.json--  watcher
 C++  --s<N>_conv.json (timings) , s<N>_solved.json (sweep)-->  (results)
```

`s<N>_reconfig.json` is left in place after being read; its presence is how the
watcher recognizes an already-answered request, and `s<N>_conv.json`'s presence
marks the request fully consumed.

### `s<N>_request.json` — C++ → watcher (the surrogate)
The reduced operator the watcher uses to choose orderings. Sized by block count
and Krylov rank only (never the global problem size `n`).
```json
{
  "n_blocks": 3,
  "block_sizes": [8, 6, 4],
  "krylov_dim": 5,                 // curDim (Krylov vectors captured)
  "ranks": [r0, r1, r2],           // per-block truncated rank, <= krylov_dim
  "equation_ends": [0, r0, r0+r1, R],  // block boundaries in C_hat/b_hat, R=sum(ranks)
  "C_hat": [[...], ...],           // R x R reduced operator
  "b_hat": [...]                   // length-R reduced RHS
}
```

### `s<N>_reconfig.json` — watcher → C++ (the response)
The watcher's answer. C++ reads `use_ordering` (required), `selection_mode`
(optional, default `"chosen"`), and `test_orderings` (optional); `opt_ordering`
and `exh_opt_ordering` are written but not consumed by C++ (reserved for
recording the surrogate-search vs exhaustive-search picks).
```json
{
  "selection_mode": "chosen",      // "chosen" | "best_conv" | "best_time"
  "use_ordering": [0, 0, 1],       // applied in "chosen" mode
  "opt_ordering": [0, 0, 1],       // record-only (surrogate's pick)
  "exh_opt_ordering": [],          // record-only (exhaustive pick; blank for now)
  "test_orderings": [[0,1,2],[0,0,1],[0,0,0]]  // candidates to build/solve/time; may be []
}
```
An **ordering** is a restricted-growth vector: `ordering[k] = g` puts original
block `k` into new group `g` (group ids contiguous from 0). A single-member
group is kept as-is; a multi-member group is assembled into one monolithic
block and factored jointly. `selection_mode` decides which ordering is applied
and returned:
- **`chosen`** — apply `use_ordering`. The `test_orderings` sweep, if present,
  is still built/timed and recorded, but `use_ordering` is what's returned.
- **`best_conv`** — sweep all `test_orderings`, return the one with the best
  convergence (fewest iterations; ties → lower residual; converged first).
- **`best_time`** — sweep all `test_orderings`, return the fastest
  (factor + iterate wall time) among those that converged.

### `s<N>_conv.json` — C++ output (per-solve timings)
The first solve vs. the applied re-solve. `solve2` is `null` if no re-solve
happened (e.g. the watcher timed out). On a rescued stall, `solve1` is the
stall (`iterations` == max iters) and `solve2` the converged rescue.
```json
{
  "request_id": 2,
  "solve1": {
    "iterations": 25,
    "initial_residual": 9.05, "final_residual": 1.1e-8,
    "factor_wall_time_sec": 0.41, "iterate_wall_time_sec": 0.55,
    "total_wall_time_sec": 0.96, "wall_time_sec": 0.96  // == total (legacy alias)
  },
  "solve2": { ... same shape ... }   // or null
}
```

### `s<N>_solved.json` — C++ output (the sweep, only when `test_orderings` given)
Every candidate ordering built/solved/timed on the as-if-original flat system,
plus which one was selected. `selected_index` is the row in `results` that was
returned (or -1 if the selected ordering — e.g. `use_ordering` in `chosen`
mode — isn't among the candidates).
```json
{
  "selection_mode": "best_conv",
  "request_id": 2,
  "selected_index": 0,
  "selected_ordering": [0, 0, 0],
  "results": [
    { "ordering": [0,0,0], "iterations": 1, "converged": true,
      "initial_residual": 9.05, "final_residual": 1.5e-16,
      "factor_wall_time_sec": 0.11, "iterate_wall_time_sec": 0.09,
      "total_wall_time_sec": 0.20 },
    ...
  ]
}
```
All wall times in `conv.json`/`solved.json` are the cross-rank maximum
(critical path), so they're identical on every rank and time-based selection
is deterministic.
