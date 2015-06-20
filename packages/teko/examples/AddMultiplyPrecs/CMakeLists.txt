
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
   AddMultiply_Example
   SOURCES 
      Driver.cpp
   COMM mpi
   )

TRIBITS_COPY_FILES_TO_BINARY_DIR(addmultiply_example_copyfiles
   SOURCE_FILES
      solverparams.xml
      lhs_test.mm
      rhs_test.mm
      modified.mm
   )
