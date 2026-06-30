IF(${PACKAGE_NAME}_ENABLE_Stratimikos
   AND ${PACKAGE_NAME}_ENABLE_Thyra
   AND ${PACKAGE_NAME}_ENABLE_Teko
   AND ${PACKAGE_NAME}_ENABLE_Belos
)
  TRIBITS_ADD_EXECUTABLE(
    Teko
    SOURCES Teko.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    InverseFactoryOperator-Tpetra
    SOURCES inverse-factory-operator.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    StratimikosSolve
    SOURCES stratimikos-solve.cpp
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    BasicExamples
    SOURCE_FILES simple_example.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    EXEDEPS InverseFactoryOperator-Tpetra
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    BasicExamples-Strat
    SOURCE_FILES strat_example.xml
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
    EXEDEPS StratimikosSolve
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    example_copyfiles
    SOURCE_FILES nslhs_test.mm nsrhs_test.mm
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../teko/examples/data/"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
    EXEDEPS InverseFactoryOperator-Tpetra StratimikosSolve
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    example_copyfiles_from_tests
    SOURCE_FILES nsjac.mm
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../teko/tests/data/"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/data"
    EXEDEPS InverseFactoryOperator-Tpetra StratimikosSolve
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    InverseFactoryOperator-Tpetra
    NAME "MueLuTekoInverseFactoryOperator-Tpetra"
    NUM_MPI_PROCS 1
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    StratimikosSolve
    NAME "MueLuTekoStratimikosSolve"
    NUM_MPI_PROCS 1
    COMM serial mpi
  )

ENDIF()
