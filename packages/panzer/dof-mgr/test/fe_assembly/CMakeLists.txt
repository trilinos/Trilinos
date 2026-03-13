TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(UNIT_TEST_DRIVER ${PANZER_UNIT_TEST_MAIN})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_FE_Assembly_HEX
  SOURCES test_fe_assembly_HEX.cpp ../cartesian_topology/CartesianConnManager.cpp
  NUM_MPI_PROCS 2
  COMM mpi
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Test_FE_Projection
  SOURCES test_fe_projection.cpp ../cartesian_topology/CartesianConnManager.cpp
  NUM_MPI_PROCS 2
  COMM mpi
  ARGS "--shape=Hexahedron --space=HGrad" "--shape=Tetrahedron --space=HGrad" "--shape=Quadrilateral --space=HGrad" "--shape=Triangle --space=HGrad"
    "--shape=Hexahedron --space=HCurl" "--shape=Tetrahedron --space=HCurl" "--shape=Quadrilateral --space=HCurl" "--shape=Triangle --space=HCurl"
    "--shape=Hexahedron --space=HDiv" "--shape=Tetrahedron --space=HDiv" "--shape=Quadrilateral --space=HDiv" "--shape=Triangle --space=HDiv"
    "--shape=Hexahedron --space=HVol" "--shape=Tetrahedron --space=HVol" "--shape=Quadrilateral --space=HVol" "--shape=Triangle --space=HVol"
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )
