INCLUDE(TribitsAddExecutableAndTest)   

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES_TEST Main.cpp TestFixture.cpp)

IF( Kokkos_ENABLE_Cuda )

  # pick up subpackage include path
  CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS} )
  CUDA_COMPILE( OBJECTS_CUDA TestFixture.cu )

  LIST( APPEND SOURCES_TEST ${OBJECTS_CUDA} )

ENDIF()

# The Example subpackage (of which this is a part) defines
# a dependency on KokkosCore so all KokkosCore libraries
# will automatically be added.
#
TRIBITS_ADD_LIBRARY(
  kokkosexample_fixture
  NOINSTALLHEADERS BoxElemFixture.hpp BoxElemPart.hpp HexElement.hpp
  SOURCES BoxElemPart.cpp
  )

# Automatically picks up 'kokkosexample_fixture'
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TestFixture
  SOURCES ${SOURCES_TEST}
  )

