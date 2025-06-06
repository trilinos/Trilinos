ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

IF (Sacado_ENABLE_Kokkos AND Sacado_ENABLE_KokkosKernels AND Sacado_ENABLE_TeuchosCore)

  TRIBITS_ADD_EXECUTABLE(
    FadAssemblyViewTest
    SOURCES BoxElemFixture.hpp
            BoxElemPart.hpp
            BoxElemPart.cpp
            HexElement.hpp
            fenl_functors.hpp
            TestAssembly.hpp
            TestAssembly.cpp
    COMM serial mpi
    )

ENDIF()
