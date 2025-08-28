TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_SinCos
  SOURCES Tempus_Test_BackwardEuler_SinCos.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

IF (TEMPUS_ENABLE_EPETRA_STACK OR TEMPUS_ENABLE_TPETRA_STACK)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BackwardEuler_CDR
    SOURCES Tempus_Test_BackwardEuler_CDR.cpp ${TEMPUS_UNIT_TEST_MAIN}
    TESTONLYLIBS tempus_test_models
    NUM_MPI_PROCS 1
    RUN_SERIAL
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_VanDerPol
  SOURCES Tempus_Test_BackwardEuler_VanDerPol.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_OptInterface
  SOURCES Tempus_Test_BackwardEuler_OptInterface.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_Combined_FSA
  SOURCES Tempus_BackwardEuler_Combined_FSA.cpp Tempus_BackwardEuler_FSA.hpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_Staggered_FSA
  SOURCES Tempus_BackwardEuler_Staggered_FSA.cpp Tempus_BackwardEuler_FSA.hpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_PseudoTransient_SA
  SOURCES Tempus_BackwardEuler_PseudoTransient_SA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  BackwardEuler_ASA
  SOURCES Tempus_BackwardEuler_ASA.cpp ${TEMPUS_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_SinCos_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SinCos.xml
  EXEDEPS BackwardEuler_SinCos
  )

IF (TEMPUS_ENABLE_EPETRA_STACK OR TEMPUS_ENABLE_TPETRA_STACK)
  TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_CDR_CopyFiles
    DEST_FILES Tempus_BackwardEuler_CDR.xml
    EXEDEPS BackwardEuler_CDR
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_VanDerPol_CopyFiles
  DEST_FILES Tempus_BackwardEuler_VanDerPol.xml
  EXEDEPS BackwardEuler_VanDerPol
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_OptInterface_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SinCos.xml
  EXEDEPS BackwardEuler_OptInterface
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_ASA_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SinCos_ASA.xml
  EXEDEPS BackwardEuler_ASA
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_BackwardEuler_PseudoTransient_SA_CopyFiles
  DEST_FILES Tempus_BackwardEuler_SteadyQuadratic.xml
  EXEDEPS BackwardEuler_PseudoTransient_SA
  )
