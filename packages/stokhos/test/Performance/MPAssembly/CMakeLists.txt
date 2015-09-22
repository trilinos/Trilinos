ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_ADD_EXECUTABLE(
  MPVectorAssemblyTest
  SOURCES BoxElemFixture.hpp
          BoxElemPart.hpp
          BoxElemPart.cpp
          HexElement.hpp
          fenl_functors.hpp
          TestAssembly.hpp
          TestAssembly.cpp
  COMM serial mpi
  )
