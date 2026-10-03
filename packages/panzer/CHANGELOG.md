# Change Log

## Unreleased (in develop branch)

- Panzer uses the unified Phalanx space names, `PHX::ExecutionSpace` and
  `PHX::MemorySpace`, and no longer relies on `PHX::Device` being an execution
  space. Intrepid2 and `Tpetra` template arguments that wanted a device were
  moved to the device position, where those templates have always wanted one.

- `panzer::TpetraNodeType` now pairs `PHX::ExecutionSpace` with
  `PHX::MemorySpace`, so Tpetra objects Panzer creates allocate where Phalanx
  fields allocate. Every node type and node template argument in Panzer goes
  through that typedef; nothing else names a node or a space in a node slot.
  The two forms are the same type unless shared space is enabled, so this
  changes nothing for existing configurations.

- A `static_assert` in `Panzer_NodeType.hpp` checks that Tpetra instantiated
  the resulting node. **`Tpetra_INST_HIP` and `Tpetra_INST_SYCL` are OFF by
  default**, while Kokkos picks HIP or SYCL as its default execution space, so
  a GPU build that does not set the matching option gets a host node from
  Tpetra and fails to link.

- `Phalanx_ENABLE_SHARED_SPACE` and `Tpetra_ALLOCATE_IN_SHARED_SPACE` must
  agree on a GPU build. Both select the same `Kokkos::SharedSpace`, so enabling
  one alone produces a node type Tpetra did not instantiate; the `static_assert`
  above catches it.

- `panzer::createIntrepid2Basis` and `panzer::PureBasis::getIntrepid2Basis`
  name their first template parameter `DeviceType`, which is what it has always
  been forwarded into.

- Design notes: `doc/design_notes/PanzerTpetraNodeTypePlan.txt`, whose decisions
  register records what was chosen and what was rejected for node consistency.

## Trilinos 17.2.1
