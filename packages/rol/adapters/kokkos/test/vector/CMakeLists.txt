
IF(ROL_ENABLE_Kokkos AND Kokkos_ENABLE_CUDA)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    KokkosVectorInterface
    SOURCES test_cuda.cpp
    ARGS PrintItAll
    COMM serial mpi  
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()


