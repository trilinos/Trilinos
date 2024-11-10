TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SIERRA_GRAD_PROJECTION_SOURCES
  SierraGradProjection.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  SierraGradProjection
  SOURCES ${SIERRA_GRAD_PROJECTION_SOURCES}
  CATEGORIES BASIC PERFORMANCE
  )

# Performance test representing gradient projection as used in Sierra 
TRIBITS_ADD_TEST(
   SierraGradProjection
   NAME "SierraGradProjection"
   ARGS ""
   COMM mpi
   NUM_MPI_PROCS 1
   CATEGORIES PERFORMANCE
   RUN_SERIAL
   )

