
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  TRIBITS_INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  ROL_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_02
    SOURCES example_02.cpp
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_03
    SOURCES example_03.cpp
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_04
    SOURCES example_04.cpp
    ADD_DIR_TO_NAME
  )

  ROL_ADD_EXECUTABLE(
    example_05
    SOURCES example_05.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_TEST(
    example_01
    ARGS PrintItAll
    NUM_MPI_PROCS 4
    NUM_TOTAL_CORES_USED 4
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
  )

  #TRIBITS_ADD_TEST(
  #  example_02
  #  ARGS PrintItAll
  #  NUM_MPI_PROCS 4
  #  NUM_TOTAL_CORES_USED 4
  #  PASS_REGULAR_EXPRESSION "TEST PASSED"
  #  ADD_DIR_TO_NAME
  #)

  ROL_COPY_FILES_TO_BINARY_DIR(
    NavierStokesDataCopy
    SOURCE_FILES
      input.xml input_ex03.xml input_ex04.xml input_ex05.xml plotresults.m plotexample02.m process.sh
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

  ROL_ADD_EXECUTABLE(
    example_01_Kokkos
    SOURCES example_01K.cpp
    ADD_DIR_TO_NAME
  )

  #  ROL_ADD_EXECUTABLE(
  #    example_02_Kokkos
  #    SOURCES example_02K.cpp
  #    ADD_DIR_TO_NAME
  #  )

  #  ROL_ADD_EXECUTABLE(
  #    example_03_Kokkos
  #    SOURCES example_03K.cpp
  #    ADD_DIR_TO_NAME
  #  )
  #
  #  ROL_ADD_EXECUTABLE(
  #    example_04_Kokkos
  #    SOURCES example_04K.cpp
  #    ADD_DIR_TO_NAME
  #  )
  #
  ROL_ADD_EXECUTABLE(
    example_05_Kokkos
    SOURCES example_05K.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_TEST(
    example_01_Kokkos
    ARGS PrintItAll
    NUM_MPI_PROCS 2
    NUM_TOTAL_CORES_USED 2
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    RUN_SERIAL
  )

  #TRIBITS_ADD_TEST(
  #  example_02
  #  ARGS PrintItAll
  #  NUM_MPI_PROCS 4
  #  NUM_TOTAL_CORES_USED 4
  #  PASS_REGULAR_EXPRESSION "TEST PASSED"
  #  ADD_DIR_TO_NAME
  #)

  ROL_COPY_FILES_TO_BINARY_DIR(
    NavierStokesDataCopyK
    SOURCE_FILES
      input.xml input_ex03.xml input_ex04.xml input_ex05.xml plotresults.m plotexample02.m process.sh
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )

ENDIF()
