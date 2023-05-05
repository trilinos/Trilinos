

KOKKOS_SUBPACKAGE(Core)

IF (NOT Kokkos_INSTALL_TESTING)
  ADD_SUBDIRECTORY(src)
ENDIF()

KOKKOS_ADD_TEST_DIRECTORIES(unit_test)
  # We are using the githash etc in here, which does not work correct in Trilinos
  KOKKOS_ADD_TEST_DIRECTORIES(perf_test)

