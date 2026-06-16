# Interface & HPC issues for pairing this with a production code

Scope: bolting a real production Trilinos application onto this research hack
(the Belos adaptive hook + Krylov surrogate + reconfig handshake), for
**research runs only**. The surrogate math and timing are reusable; the list
below is what tends to bite at the seam and at scale. Ordered by how likely it
is to actually block a run.

## Interface (the production-app ↔ hack seam)

- **Type lock-in (most likely blocker).** Everything is hardcoded to
  `SC=double`, `LO=int`, `GO=long long`,
  `Node=Tpetra::KokkosCompat::KokkosOpenMPWrapperNode`
  (`Teko_KrylovSurrogate.hpp` aliases) and gated by an `if constexpr` on
  `double` / `Thyra::MultiVectorBase<double>` / `Thyra::LinearOpBase<double>`
  in `BelosBlockGmresSolMgr`. If the production code uses a different `GO` or
  `Node`, the hook either never fires (constexpr guard) or the block
  `dynamic_cast`s return null and blocks are silently treated as zero — a wrong
  surrogate with **no error**. Expect to retemplate the surrogate on the app's
  scalar/ordinal/node types.

- **GPU / non-host nodes.** The surrogate and shim read block data with host
  `getData()` loops and raw pointers. On a CUDA/HIP node that needs explicit
  host mirrors / deep copies; as written it assumes host-accessible memory.

- **Global, always-on hook.** *(mitigated)* `registerHook` still stores a
  single process-wide callback, but the hook is now opt-in: it returns
  immediately unless `TEKO_ADAPTIVE_RECONFIG` is set, so merely linking Teko no
  longer hijacks unrelated flexible blocked solves. The gate is still
  process-wide (env), not per-solve — a per-solve `ParameterList` flag isn't
  used because Belos validates the solver parameter list and rejects unknown
  keys. Fine for research; true per-solve control would need a Belos change.

- **Operator/block structure assumptions.** Blocks must be
  `Thyra::TpetraLinearOp` over `Tpetra::CrsMatrix`, and the surrogate / merged
  assembly assume each block owns an **independent, 0-based global index
  space** (true because the demo builds each block fresh). A production blocked
  operator carved from one monolithic system may have non-zero-based or
  overlapping/ghosted column maps → mis-assembled merged blocks.

- **Diagnostics mismatch.** One `solve()` secretly runs warm-up + sweep + a
  final re-solve and overwrites the LHS, but the `SolMgr`'s `getNumIters()` /
  `achievedTol()` still report **solve 1**, not the solution actually returned.
  Any app logic that reads those will be inconsistent with the answer.

- **Hardcoded paths.** `kDefaultRequestsDir` is an absolute path baked into the
  header; the app must export `TEKO_RECONFIG_REQUESTS_DIR` or writes go nowhere
  useful. The gs/jacobi "method" is inferred by a `dynamic_cast` heuristic on
  the right-preconditioner type — a different prec structure mislabels it.

## HPC / parallel

- **Filesystem + external-process handshake.** The solve blocks (up to 600 s,
  200 ms poll) waiting for a watcher to write `s<N>_reconfig.json`. At scale
  this needs a shared filesystem visible to **rank 0 and the watcher**, which
  on many clusters means the watcher can't simply live on the launch node;
  parallel-FS metadata latency and the per-solve JSON round-trip add up if
  there are many solves. For research this is usable but slow and fragile; an
  in-process selector callback would remove it entirely.

- **Rank-divergent failure → deadlock.** *(mitigated)* `solveOrdering` now
  wraps the preconditioner build, max-reduces a failure flag, and has all ranks
  skip a candidate together before entering the collective solve — so a
  singular merged block in one swept ordering no longer hangs the run. Residual
  risk: this assumes the build throws *after* its own collectives complete (the
  same assumption `teko_ext`'s working diagonal-block guard relies on); a
  failure that hangs *inside* a collective is still unhandled.

- **Watcher liveness / placement.** *(improved)* The watcher is now started
  once by the interface and stopped at interpreter exit, with a 1000 s idle and
  parent-death self-terminate as a safety net — no more per-solve spawn churn or
  orphans. The underlying need for a shared filesystem visible to rank 0 and
  the watcher remains (inherent to the file handshake); if it isn't running,
  instrumented solves still stall the full 600 s then return the first solve.

- **Cost multiplier.** *(partly reduced)* With `test_orderings` set, one
  `solve()` is still warm-up + Bell(nb) candidate solves + a final solve, but
  the Stratimikos `InverseLibrary` is now built once per hook and reused across
  all candidates instead of rebuilt per ordering. CRS blocks are still
  reassembled per merged group. Fine for small research problems; budget for it
  before pointing it at a large operator.

- **Threading.** `FactorTimeRegistry` is a global atomic and the recursion
  guard is `thread_local`; concurrent solves on multiple host threads
  cross-contaminate the factor-time accounting. Single-threaded solve is
  assumed.

- **Per-process state across solves.** Factor-time leaks between solves unless
  the registry is reset at each solve start; the `teko_ext` shim does this, a
  drop-in app must too (`Teko::FactorTimeRegistry::reset()`).

## Lowest-friction integration path
Retemplate on the app's `Node`/`GO`; replace the filesystem/Python handshake
with an in-process ordering-selector callback; gate the hook so only the
solves you care about are instrumented; and verify the block operators present
independent per-block maps before trusting merged-group assembly.
