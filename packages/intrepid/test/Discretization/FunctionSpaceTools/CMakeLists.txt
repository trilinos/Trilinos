INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_01
  SOURCES test_01.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_02
  SOURCES test_02.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  Test_03
  SOURCES test_03.cpp
  ARGS PrintItAll
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  ADD_DIR_TO_NAME
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(FunctionSpaceToolsTestDataCopy
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

