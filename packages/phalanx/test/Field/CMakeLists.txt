

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  field
  SOURCES Field.cpp
  TESTONLYLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield
  SOURCES MDField.cpp
  TESTONLYLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_runtime
  SOURCES MDField_Runtime.cpp
  TESTONLYLIBS phalanx_test_utilities
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Test passed!"
  )

IF (NOT Phalanx_ENABLE_COMPILETIME_ARRAY_CHECK)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mdfield_simulated_intrepid
    SOURCES MDField_SimulatedIntrepid.cpp
    TESTONLYLIBS phalanx_test_utilities
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Test passed!"
    )
ENDIF()

# Will add this test in a later commit

#ASSERT_DEFINED(Phalanx_ENABLE_Intrepid)
#
#IF (Phalanx_ENABLE_Intrepid)
#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    mdfield_intrepid
#    SOURCES MDField_Intrepid.cpp
#    TESTONLYLIBS phalanx_test_utilities
#    COMM serial mpi
#    NUM_MPI_PROCS 1
#    PASS_REGULAR_EXPRESSION "Test passed!"
#    )
#
#  TRIBITS_COPY_FILES_TO_BINARY_DIR(FunctionSpaceToolsTestDataCopy
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
