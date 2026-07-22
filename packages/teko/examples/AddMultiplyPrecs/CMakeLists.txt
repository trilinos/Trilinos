TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  AddMultiply_Example
  SOURCES Driver.cpp
  COMM serial mpi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  addmultiply_example_copyfiles
  SOURCE_FILES
    lhs_test.mm
    modified.mm
    rhs_test.mm
    solverparams.xml
)
