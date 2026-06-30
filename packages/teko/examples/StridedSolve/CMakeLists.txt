
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   StridedSolve_Example
   SOURCES 
      strided_solve.cpp
   COMM mpi serial
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(strided_solve_example_copyfiles
  SOURCE_FILES solverparams.xml
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles
  SOURCE_FILES nslhs_test.mm
  nsrhs_test.mm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/../data"
  EXEDEPS StridedSolve_Example
)

IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(example_copyfiles_from_tests
    SOURCE_FILES nsjac.mm
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../tests/data"
    DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/../data"
    EXEDEPS StridedSolve_Example
  )
ENDIF()
