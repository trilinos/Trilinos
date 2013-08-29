INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(LIBRARIES kokkosexample_fixture)

# The Example subpackage (of which this is a part) defines
# a dependencie on KokkosCore so all KokkosCore libraries
# will automatically be added.
#
TRIBITS_ADD_LIBRARY(
  kokkosexample_fixture
  NOINSTALLHEADERS BoxElemFixture.hpp BoxElemPart.hpp HexElement.hpp
  SOURCES BoxElemPart.cpp
  )

IF( Kokkos_ENABLE_Cuda )

  # The 'CUDALIBRARY' tells Tribits to use the
  # cuda compiler and library archiver
  #
  TRIBITS_ADD_LIBRARY(
    kokkosexample_fixture_cuda
    NOINSTALLHEADERS BoxElemFixture.hpp BoxElemPart.hpp HexElement.hpp
    SOURCES TestFixture.cu
    CUDALIBRARY
    )

  LIST( APPEND LIBRARIES
    kokkosexample_fixture_cuda
  )

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestFixture
  SOURCES Main.cpp TestFixture.cpp TestFixture.cu
  )

