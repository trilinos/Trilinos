INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)


INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../thyra/core/example/operator_solve)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Belos)
   PACKAGE_ADD_EXECUTABLE_AND_TEST(
     Belos_GCRODR_strattest 
     SOURCES
       test_belos_gcrodr.cpp
       ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
     COMM serial mpi
     NUM_MPI_PROCS 4
     )
ENDIF()

PACKAGE_ADD_EXECUTABLE(
  test_single_stratimikos_solver_driver
  SOURCES
    test_single_stratimikos_solver_driver.cpp
    test_single_stratimikos_solver.hpp
    test_single_stratimikos_solver.cpp
  COMM serial mpi
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestFourByFourMtx
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS test_single_stratimikos_solver_driver
  )


PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestXMLFiles
  DEST_FILES
    FourByFour.xml
    FourByFour.amesos.xml
    FourByFour.aztecoo.np.xml
    FourByFour.aztecoo.ifpack.xml
    FourByFour.aztecoo.ifpack.adjoint.xml
    FourByFour.aztecoo.ml.xml
    FourByFour.belos.np.xml
    FourByFour.belos.ifpack.xml
    FourByFour.belos.ml.xml
    BelosGCRODRTest.xml
  SOURCE_PREFIX "_"
  EXEDEPS test_single_stratimikos_solver_driver
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Amesos)
  PACKAGE_ADD_TEST(
    test_single_stratimikos_solver_driver
    NAME test_single_stratimikos_solver_driver_amesos
    ARGS 
      "--input-file=FourByFour.xml"
      "--input-file=FourByFour.amesos.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
IF (${PACKAGE_NAME}_ENABLE_AztecOO)


  PACKAGE_ADD_TEST(
    test_single_stratimikos_solver_driver
    NAME test_single_stratimikos_solver_driver_aztecoo_np
    ARGS "--input-file=FourByFour.aztecoo.np.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_TEST(
      test_single_stratimikos_solver_driver
      NAME test_single_stratimikos_solver_driver_aztecoo_ifpack
      ARGS
        "--input-file=FourByFour.aztecoo.ifpack.xml"
        "--input-file=FourByFour.aztecoo.ifpack.adjoint.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
  IF (${PACKAGE_NAME}_ENABLE_ML)

    PACKAGE_ADD_TEST(
      test_single_stratimikos_solver_driver
      NAME test_single_stratimikos_solver_driver_aztecoo_ml
      ARGS
        "--input-file=FourByFour.aztecoo.ml.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
IF (${PACKAGE_NAME}_ENABLE_Belos)


  PACKAGE_ADD_TEST(
    test_single_stratimikos_solver_driver
    NAME test_single_stratimikos_solver_driver_belos_np
    ARGS "--input-file=FourByFour.belos.np.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_TEST(
      test_single_stratimikos_solver_driver
      NAME test_single_stratimikos_solver_driver_belos_ifpack
      ARGS
        "--input-file=FourByFour.belos.ifpack.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
  IF (${PACKAGE_NAME}_ENABLE_ML)

    PACKAGE_ADD_TEST(
      test_single_stratimikos_solver_driver
      NAME test_single_stratimikos_solver_driver_belos_ml
      ARGS
        "--input-file=FourByFour.belos.ml.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


ENDIF()
