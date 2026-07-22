TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# value types
SET(Intrepid2_TEST_ETI_VALUETYPE_NAME "")
SET(Intrepid2_TEST_ETI_VALUETYPE "")

LIST(APPEND Intrepid2_TEST_ETI_VALUETYPE_NAME "DOUBLE")
LIST(APPEND Intrepid2_TEST_ETI_VALUETYPE      "double")
LIST(APPEND Intrepid2_TEST_ETI_SACADO         "0")

LIST(LENGTH Intrepid2_TEST_ETI_VALUETYPE_NAME ETI_VALUETYPE_COUNT)
MATH(EXPR ETI_VALUETYPE_COUNT "${ETI_VALUETYPE_COUNT}-1")

# Host test
SET(Intrepid2_TEST_ETI_FILE "")
LIST(APPEND Intrepid2_TEST_ETI_FILE 
  "test_orientation_encoding"
  "test_orientation_HEX"
  "test_orientation_HEX_newBasis"
  "test_orientation_PYR_QuadFace_newBasis"
  "test_orientation_PYR_TriFace_newBasis"
  "test_orientation_QUAD"
  "test_orientation_QUAD_newBasis"
  "test_orientation_TET"
  "test_orientationtools_coeff_matrix"
  "test_orientationtools_mesh"
  "test_orientationtools_modify_basis_hcurl"
  "test_orientationtools_modify_basis_hdiv"
  "test_orientationtools_modify_basis_hgrad"
  "test_orientationtools_modify_basis_lower_order"
  "test_orientation_TRI")

SET(Intrepid2_TEST_ETI_DEVICE_NAME "")
SET(Intrepid2_TEST_ETI_DEVICE "")
IF(Kokkos_ENABLE_SERIAL)
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "Serial")
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>")
ENDIF()

#  disable openmp for PR test failure due to wall time
#IF(Kokkos_ENABLE_OPENMP)
#  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "OpenMP")
#  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>")
#ENDIF()

LIST(LENGTH Intrepid2_TEST_ETI_DEVICE_NAME ETI_DEVICE_COUNT)
MATH(EXPR ETI_DEVICE_COUNT "${ETI_DEVICE_COUNT}-1")

IF (${ETI_DEVICE_COUNT} GREATER_EQUAL 0)
  FOREACH(I RANGE ${ETI_DEVICE_COUNT})
    LIST(GET Intrepid2_TEST_ETI_DEVICE_NAME ${I} ETI_DEVICE_NAME)
    LIST(GET Intrepid2_TEST_ETI_DEVICE      ${I} ETI_DEVICE)
    FOREACH(J RANGE ${ETI_VALUETYPE_COUNT})
      LIST(GET Intrepid2_TEST_ETI_VALUETYPE_NAME ${J} ETI_VALUETYPE_NAME)
      LIST(GET Intrepid2_TEST_ETI_VALUETYPE      ${J} ETI_VALUETYPE)
      LIST(GET Intrepid2_TEST_ETI_SACADO         ${J} ETI_SACADO)
      FOREACH(ETI_FILE IN LISTS Intrepid2_TEST_ETI_FILE)
        SET(ETI_NAME "${ETI_FILE}_${ETI_DEVICE_NAME}_${ETI_VALUETYPE_NAME}")
        CONFIGURE_FILE(eti/${ETI_FILE}_ETI.in ${ETI_NAME}.cpp)
  
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          ${ETI_NAME}
          SOURCES ${ETI_NAME}.cpp 
          ARGS PrintItAll
          NUM_MPI_PROCS 1
          PASS_REGULAR_EXPRESSION "TEST PASSED"
          ADD_DIR_TO_NAME
        )
  
      ENDFOREACH()
    ENDFOREACH()
  ENDFOREACH()
ENDIF()

# Device test
SET(Intrepid2_TEST_ETI_FILE "")
LIST(APPEND Intrepid2_TEST_ETI_FILE
  "test_orientationtools_modify_basis_hcurl"
  "test_orientationtools_modify_basis_hdiv"
  "test_orientationtools_modify_basis_hgrad"
  "test_orientationtools_modify_basis_lower_order")

SET(Intrepid2_TEST_ETI_DEVICE_NAME "")
SET(Intrepid2_TEST_ETI_DEVICE "")

IF(Kokkos_ENABLE_CUDA)
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "CUDA")
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>")
ENDIF()
IF(Kokkos_ENABLE_HIP)
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "HIP")
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::HIP,Kokkos::HIPSpace>")
ENDIF()

LIST(LENGTH Intrepid2_TEST_ETI_DEVICE_NAME ETI_DEVICE_COUNT)
MATH(EXPR ETI_DEVICE_COUNT "${ETI_DEVICE_COUNT}-1")

IF (${ETI_DEVICE_COUNT} GREATER_EQUAL 0)
  FOREACH(I RANGE ${ETI_DEVICE_COUNT})
    LIST(GET Intrepid2_TEST_ETI_DEVICE_NAME ${I} ETI_DEVICE_NAME)
    LIST(GET Intrepid2_TEST_ETI_DEVICE      ${I} ETI_DEVICE)
    FOREACH(J RANGE ${ETI_VALUETYPE_COUNT})
      LIST(GET Intrepid2_TEST_ETI_VALUETYPE_NAME ${J} ETI_VALUETYPE_NAME)
      LIST(GET Intrepid2_TEST_ETI_VALUETYPE      ${J} ETI_VALUETYPE)
      LIST(GET Intrepid2_TEST_ETI_SACADO         ${J} ETI_SACADO)
      FOREACH(ETI_FILE IN LISTS Intrepid2_TEST_ETI_FILE)
        SET(ETI_NAME "${ETI_FILE}_${ETI_DEVICE_NAME}_${ETI_VALUETYPE_NAME}")
        CONFIGURE_FILE(eti/${ETI_FILE}_ETI.in ${ETI_NAME}.cpp)
  
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          ${ETI_NAME}
          SOURCES ${ETI_NAME}.cpp 
          ARGS PrintItAll
          NUM_MPI_PROCS 1
          PASS_REGULAR_EXPRESSION "TEST PASSED"
          ADD_DIR_TO_NAME
        )
  
      ENDFOREACH()
    ENDFOREACH()
  ENDFOREACH()
ENDIF()
