ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Sacado_ENABLE_KokkosCore AND Sacado_ENABLE_KokkosContainers AND
    Sacado_ENABLE_Teuchos)

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
