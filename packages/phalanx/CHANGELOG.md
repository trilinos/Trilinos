# Change Log

## Unreleased (in develop branch)

- `PHX::Device` is now `Kokkos::Device<PHX::ExecutionSpace, PHX::MemorySpace>`
  rather than a bare execution space, so that a configured memory space reaches
  the Views. A `Kokkos::Device` has no `size_type` or `array_layout` and is not
  an execution space instance, so uses of `PHX::Device` that meant an execution
  space must say `PHX::ExecutionSpace`.
  `Phalanx_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE=ON` restores the old
  meaning for applications that have not migrated.

- The execution and memory spaces are spelled `PHX::ExecutionSpace` and
  `PHX::MemorySpace`. `PHX::exec_space`, `PHX::ExecSpace`, `PHX::mem_space` and
  `PHX::MemSpace` remain as deprecated aliases. Use
  `Phalanx_SHOW_DEPRECATED_WARNINGS=OFF` to silence the warnings while
  migrating, and `Phalanx_HIDE_DEPRECATED_CODE=ON` to prove a migration is
  complete.

- `Phalanx_ENABLE_SHARED_SPACE=ON` is supported: the
  memory space becomes `Kokkos::SharedSpace`. A `static_assert` checks that the
  configured execution space can reach the memory space it is paired with.

- Removed `Phalanx_DEFAULT_MEMORY_SPACE`. The memory space is derived, never
  chosen: the execution space's own, or `Kokkos::SharedSpace`. Nothing read
  this variable before `PHX::Device` became a `Kokkos::Device`.

- Added `scripts/migrate_phx_device.py` migrates application code to the new
  device, memory space and execution space definitions, and reports what it
  cannot decide rather than guessing. Run it with `--check` first.
  `scripts/test_migrate_phx_device.py` is its test suite.

- Design notes: `doc/design_notes/PhalanxSharedSpacePlan.txt`, whose decisions
  register records what was chosen and what was rejected.

## Trilinos 17.2.1
