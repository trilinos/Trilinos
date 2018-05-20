
IF(${PROJECT_NAME}_ENABLE_Intrepid AND
   ${PROJECT_NAME}_ENABLE_Ifpack2  AND
   ${PROJECT_NAME}_ENABLE_MueLu    AND
   ${PROJECT_NAME}_ENABLE_Amesos2  AND
   ${PROJECT_NAME}_ENABLE_Tpetra )

  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/teuchos/src/sol)

  # Need ROL_TpetraMultiVector.hpp
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/vector)
  INCLUDE_DIRECTORIES(${${PACKAGE_NAME}_SOURCE_DIR}/adapters/tpetra/src/sol)

  TRIBITS_ADD_EXECUTABLE(
    example_01
    SOURCES example_01.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_02
    SOURCES example_02.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_03
    SOURCES example_03.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_04
    SOURCES example_04.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_05
    SOURCES example_05.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_06
    SOURCES example_06.cpp
    ADD_DIR_TO_NAME
  )

  TRIBITS_ADD_EXECUTABLE(
    example_07
    SOURCES example_07.cpp
    ADD_DIR_TO_NAME
  )

#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    example_02
#    SOURCES example_02.cpp
#    ARGS PrintItAll
#    NUM_MPI_PROCS 4
#    NUM_TOTAL_CORES_USED 4
#    PASS_REGULAR_EXPRESSION "TEST PASSED"
#    ADD_DIR_TO_NAME
#  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    ThermalFluidsDataCopy
    SOURCE_FILES
      input.xml input_ex03.xml input_ex04.xml input_ex05.xml input_ex06.xml input_ex07.xml tf-cube-3x3x3.txt tf-cube-15x15x15.txt plotresults.m plotresults_ex03.m plotresults_ex04.m process3D.sh
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  )


ENDIF()
