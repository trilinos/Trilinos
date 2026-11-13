ADD_SUBDIRECTORY(VerySimple)

IF(${PACKAGE_NAME}_ENABLE_Xpetra AND Tpetra_INST_INT_LONG_LONG)
  TRIBITS_ADD_EXECUTABLE(
    CreateCrsGraphFromStencil
    SOURCES CreateCrsGraphFromStencil.cpp
    COMM serial mpi
  )

  FOREACH(TEST_NAME
      Serial1DThreePoint
      Serial2DFivePoint
      Serial3DSevenPoint
      Serial2DNinePointTwoDofs)
    TRIBITS_ADD_TEST(
      CreateCrsGraphFromStencil
      NAME CreateCrsGraphFromStencil_${TEST_NAME}
      ARGS ${TEST_NAME}
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDFOREACH()

  FOREACH(TEST_NAME
      Mpi2DFivePoint
      Mpi3D27PointThreeDofs
      Mpi2DBigStarRadiusTwo)
    TRIBITS_ADD_TEST(
      CreateCrsGraphFromStencil
      NAME CreateCrsGraphFromStencil_${TEST_NAME}
      ARGS ${TEST_NAME}
      COMM mpi
      NUM_MPI_PROCS 4
      STANDARD_PASS_OUTPUT
    )
  ENDFOREACH()
ENDIF()


IF(${PACKAGE_NAME}_ENABLE_Xpetra AND NOT ${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bug_5851
      SOURCES bug_5851.cpp
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

ENDIF()

