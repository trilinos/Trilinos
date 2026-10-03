TRIBITS_SUBPACKAGE(DofMgr)

# panzer::TpetraNodeType pairs PHX::ExecutionSpace with PHX::MemorySpace, so
# Phalanx and Tpetra have to agree about shared memory.  Both options select
# the same Kokkos::SharedSpace -- CudaUVMSpace, HIPManagedSpace or
# SYCLSharedUSMSpace -- so enabling one alone gives panzer a node type Tpetra
# did not instantiate, and every Tpetra object panzer creates fails to link.
#
# Only checked when a GPU backend is enabled, and that is not timidity: with
# host backends alone Kokkos defines SharedSpace AS HostSpace, which is already
# the host execution space's own memory space, so the two nodes are the same
# type whatever the options say and there is nothing to disagree about.
#
# The authority is the static_assert in Panzer_NodeType.hpp, which compares the
# actual types.  This check exists only to say so at configure time instead of
# part way through a build.  Deliberately NOT checked here: whether Tpetra
# instantiated the node at all.  That would mean restating Tpetra's node
# selection rules in panzer, where they would quietly rot the next time Tpetra
# changes them.
IF ((Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL) AND
    (NOT Phalanx_ENABLE_SHARED_SPACE STREQUAL Tpetra_ALLOCATE_IN_SHARED_SPACE))
  MESSAGE(FATAL_ERROR "ERROR (PanzerDofMgr): Phalanx_ENABLE_SHARED_SPACE=${Phalanx_ENABLE_SHARED_SPACE} and Tpetra_ALLOCATE_IN_SHARED_SPACE=${Tpetra_ALLOCATE_IN_SHARED_SPACE} disagree. On a GPU build both select Kokkos::SharedSpace, so panzer::TpetraNodeType would name a node Tpetra did not instantiate and panzer's Tpetra objects would fail to link. Set both ON or both OFF.")
ENDIF()

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()

# TOUCH (to trigger reconfigure)
