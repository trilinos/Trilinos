
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  field
  SOURCES Field.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  mdfield
  SOURCES MDField.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  mdfield_runtime
  SOURCES MDField_Runtime.cpp
  DEPLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

# Will add this test in a later commit

#ASSERT_DEFINED(Phalanx_ENABLE_Intrepid)
#
#IF (Phalanx_ENABLE_Intrepid)
#  PACKAGE_ADD_EXECUTABLE_AND_TEST(
#    mdfield_intrepid
#    SOURCES MDField_Intrepid.cpp
#    DEPLIBS phalanx_test_utilities
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    PASS_REGULAR_EXPRESSION "Test passed!"
#    )
#
#  PACKAGE_COPY_FILES_TO_BINARY_DIR(FunctionSpaceToolsTestDataCopy
#    SOURCE_FILES
#      mass_HCURL_HEX_I1_FEM_01.dat
#      mass_HCURL_HEX_I1_FEM_02.dat
#      mass_HCURL_HEX_I1_FEM_03.dat
#      mass_HCURL_HEX_I1_FEM_04.dat
#      mass_TET_FEM_P1_01.dat
#      mass_TET_FEM_P1_02.dat
#      mass_TET_FEM_P1_03.dat
#      mass_TET_FEM_P1_04.dat
#      mass_fp_HCURL_HEX_I1_FEM_01.dat
#      mass_fp_HCURL_HEX_I1_FEM_02.dat
#      mass_fp_HCURL_HEX_I1_FEM_03.dat
#      mass_fp_HCURL_HEX_I1_FEM_04.dat
#      stiff_HCURL_HEX_I1_FEM_01.dat
#    stiff_HCURL_HEX_I1_FEM_02.dat
#    stiff_HCURL_HEX_I1_FEM_03.dat
#    stiff_HCURL_HEX_I1_FEM_04.dat
#    stiff_TET_FEM_P1_01.dat
#    stiff_TET_FEM_P1_02.dat
#    stiff_TET_FEM_P1_03.dat
#    stiff_TET_FEM_P1_04.dat
#    stiff_fp_HCURL_HEX_I1_FEM_01.dat
#    stiff_fp_HCURL_HEX_I1_FEM_02.dat
#    stiff_fp_HCURL_HEX_I1_FEM_03.dat
#    stiff_fp_HCURL_HEX_I1_FEM_04.dat
#  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/testdata"
#  DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/testdata"
#)
#
#ENDIF(Phalanx_ENABLE_Intrepid)