
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_01
    SOURCES example_01.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_02
    SOURCES example_02.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_03
    SOURCES example_03.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_04
    SOURCES example_04.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    StochasticStefanBoltzmannDataCopy
    SOURCE_FILES
      input.xml
      input_app.xml
      input_ex04.xml
      SROMinput.xml
      plotresults.m
      plotstats.m
      process.sh
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()

IF(${PROJECT_NAME}_ENABLE_Intrepid2 AND
   ${PROJECT_NAME}_ENABLE_Ifpack2   AND
   ${PROJECT_NAME}_ENABLE_MueLu     AND
   ${PROJECT_NAME}_ENABLE_Amesos2   AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_01_Kokkos
    SOURCES example_01K.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_02_Kokkos
    SOURCES example_02K.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE_AND_TEST(
    example_03_Kokkos
    SOURCES example_03K.cpp
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_04_Kokkos
    SOURCES example_04K.cpp
    ADD_DIR_TO_NAME
  )

  ROL_COPY_FILES_TO_BINARY_DIR(
    StochasticStefanBoltzmannDataCopyK
    SOURCE_FILES
      input.xml
      input_app.xml
      input_ex04.xml
      SROMinput.xml
      plotresults.m
      plotstats.m
      process.sh
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
