
INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(example_SOURCES
  Example_MultiWorkset.cpp
  AlgebraicTypes.hpp
  Traits.hpp
  Traits.cpp
  MultiWorkset_FactoryTraits.hpp
  Cell.hpp
  Cell.cpp
  Edge.hpp
  MultiWorkset.hpp
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  example_multiworkset
  SOURCES ${example_SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Run has completed successfully!"
  )
