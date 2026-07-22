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
SET(Intrepid2_TEST_ETI_FILE "test_00;test_01;test_02;test_03;test_04")

SET(Intrepid2_TEST_ETI_DEVICE_NAME "")
SET(Intrepid2_TEST_ETI_DEVICE "")
IF(Kokkos_ENABLE_SERIAL)
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "Serial")
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::Serial,Kokkos::HostSpace>")
ENDIF()
IF(Kokkos_ENABLE_OPENMP)
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE_NAME "OpenMP")
  LIST(APPEND Intrepid2_TEST_ETI_DEVICE      "Kokkos::Device<Kokkos::OpenMP,Kokkos::HostSpace>")
ENDIF()
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

#
# move data files
#

TRIBITS_COPY_FILES_TO_BINARY_DIR(FunctionSpaceToolsTestDataCopy
  SOURCE_FILES
    mass_HCURL_HEX_I1_FEM_01.dat
    mass_HCURL_HEX_I1_FEM_02.dat
    mass_HCURL_HEX_I1_FEM_03.dat
    mass_HCURL_HEX_I1_FEM_04.dat
    mass_HDIV_HEX_I1_FEM_01.dat
    mass_HDIV_HEX_I1_FEM_02.dat
    mass_HDIV_HEX_I1_FEM_03.dat
    mass_HDIV_HEX_I1_FEM_04.dat
    mass_TET_FEM_P1_01.dat
    mass_TET_FEM_P1_02.dat
    mass_TET_FEM_P1_03.dat
    mass_TET_FEM_P1_04.dat
    mass_fp_HCURL_HEX_I1_FEM_01.dat
    mass_fp_HCURL_HEX_I1_FEM_02.dat
    mass_fp_HCURL_HEX_I1_FEM_03.dat
    mass_fp_HCURL_HEX_I1_FEM_04.dat
    mass_fp_HDIV_HEX_I1_FEM_01.dat
    mass_fp_HDIV_HEX_I1_FEM_02.dat
    mass_fp_HDIV_HEX_I1_FEM_03.dat
    mass_fp_HDIV_HEX_I1_FEM_04.dat
    stiff_HCURL_HEX_I1_FEM_01.dat
    stiff_HCURL_HEX_I1_FEM_02.dat
    stiff_HCURL_HEX_I1_FEM_03.dat
    stiff_HCURL_HEX_I1_FEM_04.dat
    stiff_HDIV_HEX_I1_FEM_01.dat
    stiff_HDIV_HEX_I1_FEM_02.dat
    stiff_HDIV_HEX_I1_FEM_03.dat
    stiff_HDIV_HEX_I1_FEM_04.dat
    stiff_TET_FEM_P1_01.dat
    stiff_TET_FEM_P1_02.dat
    stiff_TET_FEM_P1_03.dat
    stiff_TET_FEM_P1_04.dat
    stiff_fp_HCURL_HEX_I1_FEM_01.dat
    stiff_fp_HCURL_HEX_I1_FEM_02.dat
    stiff_fp_HCURL_HEX_I1_FEM_03.dat
    stiff_fp_HCURL_HEX_I1_FEM_04.dat
    stiff_fp_HDIV_HEX_I1_FEM_01.dat
    stiff_fp_HDIV_HEX_I1_FEM_02.dat
    stiff_fp_HDIV_HEX_I1_FEM_03.dat
    stiff_fp_HDIV_HEX_I1_FEM_04.dat
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/testdata"
  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/testdata"
#  EXEDEPS Test_01 #Test_02
)

