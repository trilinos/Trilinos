
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   InverseFactoryOperator-Tpetra
   SOURCES 
      inverse-factory-operator.cpp
   COMM mpi
   )

TRIBITS_ADD_EXECUTABLE(
   StratimikosSolve
   SOURCES 
      stratimikos-solve.cpp
   COMM mpi
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(BasicExamples
  SOURCE_FILES 
    simple_example.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS InverseFactoryOperator-Tpetra
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(BasicExamples-Strat
  SOURCE_FILES 
    strat_example.xml
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}"
  EXEDEPS StratimikosSolve
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles
  SOURCE_FILES nslhs_test.mm
  nsrhs_test.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/../data"
  EXEDEPS InverseFactoryOperator-Tpetra StratimikosSolve
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles_from_tests
  SOURCE_FILES nsjac.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../tests/data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/../data"
  EXEDEPS InverseFactoryOperator-Tpetra StratimikosSolve
)

