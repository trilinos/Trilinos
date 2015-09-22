INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_compiletime
  SOURCES MDField_Compiletime.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_runtime
  SOURCES MDField_Runtime.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_utilities
  SOURCES MDField_Utilities.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  mdfield_simulated_intrepid
  SOURCES MDField_SimulatedIntrepid.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS phalanx_test_utilities
  NUM_MPI_PROCS 1
  )

IF(Phalanx_ENABLE_Intrepid)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    mdfield_intrepid
    SOURCES MDField_Intrepid.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
    TESTONLYLIBS phalanx_test_utilities
    NUM_MPI_PROCS 1
    )
ENDIF(Phalanx_ENABLE_Intrepid)


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
