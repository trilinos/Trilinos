ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (Sacado_ENABLE_Kokkos AND Sacado_ENABLE_TeuchosCore)

  IF(NOT ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8")))

    TRIBITS_ADD_EXECUTABLE(
      FadAdvectionConstBasisHierarchicalTest
      SOURCES common.hpp
              advection.hpp advection.cpp
              advection_hierarchical.hpp advection_hierarchical.cpp
              advection_hierarchical_dfad.hpp advection_hierarchical_dfad.cpp
              driver.cpp
      COMM serial mpi
      )

  ENDIF()

ENDIF()
