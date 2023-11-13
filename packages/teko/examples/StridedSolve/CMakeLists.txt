
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
