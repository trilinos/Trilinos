
IF(ROL_ENABLE_Kokkos AND Kokkos_ENABLE_CUDA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    KokkosVector_UnaryFunctions
    SOURCES test_unaryfunctions_cuda.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    KokkosVector_BinaryFunctions
    SOURCES test_binaryfunctions_cuda.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
    KokkosVector_ReductionOps
    SOURCES test_reductionops_cuda.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()


