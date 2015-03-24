

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Example_01
  SOURCES example_01.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
#  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

IF (${PACKAGE_NAME}_ENABLE_KokkosCore)
  # mfh 24 Mar 2015: These tests use Kokkos, but Kokkos(Core) is an
  # optional dependency of Intrepid.  Thus, we must protect building
  # these tests.

  #TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #  Example_02
  #  SOURCES example_02.cpp
  #  ARGS PrintItAll
  #  NUM_MPI_PROCS 1
  ##  PASS_REGULAR_EXPRESSION "TEST PASSED"
  #  ADD_DIR_TO_NAME
  #  )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_03
    SOURCES example_03.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Example_04
    SOURCES example_04.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF ()
