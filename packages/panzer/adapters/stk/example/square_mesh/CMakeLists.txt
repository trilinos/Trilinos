
INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
SET(square_mesh_SOURCES
  square_mesh.cpp
  )

SET(square_mesh_bc_SOURCES
  square_mesh_bc.cpp
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  square_mesh
  SOURCES ${square_mesh_SOURCES}
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  square_mesh_bc
  SOURCES ${square_mesh_bc_SOURCES}
  COMM serial mpi
  )
