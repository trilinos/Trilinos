# KokkosCore_config.h
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/kokkos/core/src)
# TpetraKernels_config.h
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/tpetra/kernels/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/tpetra/core/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/tpetra/classic/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/tpetra/classic/NodeAPI)
# TeuchosCore_config.h
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/teuchos/core/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/teuchos/kokkoscompat/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/packages/teuchos/kokkoscomm/src)

# Other Kokkos includes needed
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/kokkos/core/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/kokkos/containers/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/kokkos/algorithms/src)
# Other Tpetra includes needed
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/kernels/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/kernels/src/impl)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/core/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/core/src/kokkos_refactor)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/classic/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/classic/NodeAPI)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/tpetra/classic/LinAlg)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/core/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/parameterlist/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/kokkoscompat/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/comm/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/kokkoscomm/src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/teuchos/numerics/src)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/packages/moertel/morkon/src)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(INCLUDES
    mrk_api_classes.hpp
    mrk_compute_normals.hpp
    mrk_data_types.hpp
	mrk_interface_host_side_adapter.hpp
    mrk_interface_impl.hpp
	mrk_int_types.hpp
	mrk_manager_impl.hpp
	mrk_manager_tester.hpp
    mrk_search_for_pallet_generating_faces.hpp
    morkon_unit_test_utils.hpp
    morkon_unit_test_utils_TRI3.hpp
    mrk_planar_geom.hpp
   )


##
## This is just for now.  See the contents kokkos/core/unit_test to see what should be
## done to handle the different types of devices.
##

SET(CPPSOURCES
    morkon_unit_test_utils_TRI3.cpp
    morkon_core_unit_test.cpp
    morkon_planar_geom_unit_test.cpp
   )


#TRIBITS_ADD_LIBRARY(morkon HEADERS ${INCLUDES})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  morkon_core_unit_test
  NAME morkon_core_unit_test
  SOURCES morkon_core_unit_test.cpp morkon_unit_test_utils_TRI3.cpp morkon_planar_geom_unit_test.cpp
  ARGS
  COMM mpi
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

#SET_PROPERTY(TARGET morkon_exp_exp_lib APPEND PROPERTY LINKER_LANGUAGE CXX)

