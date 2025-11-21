
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_01
    SOURCES example_01.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    PoissonDataCopy
    SOURCE_FILES
      input.xml plotresults.m p-cube-4x4x4.txt p-cube-4x4x4.e p-cube-8x8x8.txt p-cube-8x8x8.e
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()

IF(${PROJECT_NAME}_ENABLE_Intrepid2 AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    example_01_Kokkos
    SOURCES example_01_K.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    PoissonDataCopyK
    SOURCE_FILES
      input.xml plotresults.m p-cube-4x4x4.txt p-cube-4x4x4.e p-cube-8x8x8.txt p-cube-8x8x8.e
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
