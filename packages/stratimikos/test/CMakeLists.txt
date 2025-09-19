

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../thyra/core/example/operator_solve)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
   TRIBITS_ADD_EXECUTABLE_AND_TEST(
     Belos_epetra_GCRODR_strattest
     SOURCES
       test_belos_epetra_gcrodr.cpp
       ${TEUCHOS_STD_UNIT_TEST_MAIN}
     ARGS " --show-timer-summary"
     COMM serial mpi
     NUM_MPI_PROCS 4
     )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestXMLFiles
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
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestFourByFourMtx
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters)
  TRIBITS_ADD_EXECUTABLE(
    test_epetra_stratimikos_solver_driver
    SOURCES
    test_epetra_stratimikos_solver_driver.cpp
    test_epetra_stratimikos_solver.hpp
    test_epetra_stratimikos_solver.cpp
    COMM serial mpi
  )
ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Amesos)
  TRIBITS_ADD_TEST(
    test_epetra_stratimikos_solver_driver
    NAME test_epetra_stratimikos_solver_driver_amesos
    ARGS
      "--input-file=FourByFour.xml --show-timer-summary"
      "--input-file=FourByFour.amesos.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters AND ${PACKAGE_NAME}_ENABLE_AztecOO)


  TRIBITS_ADD_TEST(
    test_epetra_stratimikos_solver_driver
    NAME test_epetra_stratimikos_solver_driver_aztecoo_np
    ARGS "--input-file=FourByFour.aztecoo.np.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_TEST(
      test_epetra_stratimikos_solver_driver
      NAME test_epetra_stratimikos_solver_driver_aztecoo_ifpack
      ARGS
        "--input-file=FourByFour.aztecoo.ifpack.xml --show-timer-summary"
        "--input-file=FourByFour.aztecoo.ifpack.adjoint.xml --show-timer-summary"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
  IF (${PACKAGE_NAME}_ENABLE_ML)

    TRIBITS_ADD_TEST(
      test_epetra_stratimikos_solver_driver
      NAME test_epetra_stratimikos_solver_driver_aztecoo_ml
      ARGS
        "--input-file=FourByFour.aztecoo.ml.xml --show-timer-summary"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
IF (${PACKAGE_NAME}_ENABLE_ThyraEpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Belos)


  TRIBITS_ADD_TEST(
    test_epetra_stratimikos_solver_driver
    NAME test_epetra_stratimikos_solver_driver_belos_np
    ARGS "--input-file=FourByFour.belos.np.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_TEST(
      test_epetra_stratimikos_solver_driver
      NAME test_epetra_stratimikos_solver_driver_belos_ifpack
      ARGS
        "--input-file=FourByFour.belos.ifpack.xml --show-timer-summary"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
  IF (${PACKAGE_NAME}_ENABLE_ML)

    TRIBITS_ADD_TEST(
      test_epetra_stratimikos_solver_driver
      NAME test_epetra_stratimikos_solver_driver_belos_ml
      ARGS
        "--input-file=FourByFour.belos.ml.xml --show-timer-summary"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


ENDIF()

ENDIF()

# Stratimikos now has optional test dependencies on Ifpack2 (and thus
# indirectly, Tpetra) and ThyraTpetraAdapters.

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack2)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Tpetra)
IF (${PROJECT_NAME}_ENABLE_Tpetra)
  ASSERT_DEFINED(Tpetra_INST_DOUBLE)
  ASSERT_DEFINED(Tpetra_INST_FLOAT)
  ASSERT_DEFINED(Tpetra_INST_INT_INT)
ENDIF()

IF (${PROJECT_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)

  TRIBITS_ADD_EXECUTABLE(
    test_tpetra_stratimikos_solver_driver
    SOURCES
      test_tpetra_stratimikos_solver_driver.cpp
      test_tpetra_stratimikos_solver.hpp
      test_tpetra_stratimikos_solver.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    test_tpetra_stratimikos_solver_driver
    NAME test_tpetra_stratimikos_solver_driver
    ARGS
      "--input-file=FourByFour.belos.np.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE(
    test_tpetra_belos
    SOURCES test_tpetra_belos.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_TEST(
    test_tpetra_belos
    NAME test_tpetra_belos
    ARGS
      "--input-file=TpetraBelosTest.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    test_tpetra_belos
    NAME test_tpetra_belos_sing_reduce
    ARGS
      "--input-file=TpetraSingReduceTest.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    test_tpetra_belos
    NAME test_tpetra_belos_prec
    ARGS
      "--input-file=stratimikos_BelosPrec_ParameterList.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    test_tpetra_belos
    NAME test_tpetra_belos_half_prec
    ARGS
      "--input-file=stratimikos_HalfBelosPrec_ParameterList.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 4
    EXCLUDE_IF_NOT_TRUE ${PROJECT_NAME}_ENABLE_FLOAT
    )

  TRIBITS_ADD_TEST(
    test_tpetra_belos
    NAME test_tpetra_belos_gmres_ir
    ARGS
      "--input-file=stratimikos_GMRES_IR_ParameterList.xml --show-timer-summary"
    COMM serial mpi
    NUM_MPI_PROCS 4
    EXCLUDE_IF_NOT_TRUE ${PROJECT_NAME}_ENABLE_FLOAT
    )

  ASSERT_DEFINED(Belos_SOURCE_DIR)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTpetraTestMtxFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES bcsstk12.mtx
    EXEDEPS test_tpetra_belos
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTpetraTestXMLFiles
    DEST_FILES
      TpetraBelosTest.xml
      TpetraSingReduceTest.xml
      stratimikos_BelosPrec_ParameterList.xml
      stratimikos_HalfBelosPrec_ParameterList.xml
      stratimikos_GMRES_IR_ParameterList.xml
    EXEDEPS test_tpetra_belos
  )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND Galeri_ENABLE_Xpetra)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Stratimikos_cp
    SOURCE_FILES stratimikos_jacobi.xml stratimikos_jacobi_half.xml stratimikos_jacobi_tpetra.xml
    )

  TRIBITS_ADD_EXECUTABLE(
    GaleriXpetraDriver
    SOURCES
    galeri_xpetra_driver.cpp
    COMM serial mpi
    )

  IF (Tpetra_INST_DOUBLE)

    TRIBITS_ADD_TEST(
      GaleriXpetraDriver
      NAME "Galeri_xpetra_double_Jacobi"
      ARGS "--scalarType=double --xml=stratimikos_jacobi.xml"
      NUM_MPI_PROCS 4
      )

    TRIBITS_ADD_TEST(
      GaleriXpetraDriver
      NAME "Galeri_xpetra_double_Jacobi_Tpetra"
      ARGS "--scalarType=double --xml=stratimikos_jacobi_tpetra.xml"
      NUM_MPI_PROCS 4
      )

    IF (Tpetra_INST_FLOAT
      AND (NOT Tpetra_INST_COMPLEX_DOUBLE
	OR (Tpetra_INST_COMPLEX_DOUBLE AND Tpetra_INST_COMPLEX_FLOAT)
	)
      )
      TRIBITS_ADD_TEST(
        GaleriXpetraDriver
        NAME "Galeri_xpetra_double_Jacobi_HalfPrecision"
        ARGS "--scalarType=double --xml=stratimikos_jacobi_half.xml"
        NUM_MPI_PROCS 4
        )
    ENDIF()

  ENDIF()

  IF (Tpetra_INST_FLOAT)

    TRIBITS_ADD_TEST(
      GaleriXpetraDriver
      NAME "Galeri_xpetra_float_Jacobi"
      ARGS "--scalarType=float --xml=stratimikos_jacobi.xml"
      NUM_MPI_PROCS 4
      )

  ENDIF()

  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF (Tpetra_INST_COMPLEX_DOUBLE AND Teuchos_ENABLE_COMPLEX)

    TRIBITS_ADD_TEST(
      GaleriXpetraDriver
      NAME "Galeri_xpetra_complex_double_Jacobi"
      ARGS "--scalarType=\"complex<double>\" --xml=stratimikos_jacobi.xml"
      NUM_MPI_PROCS 4
      )

  ENDIF()

  IF (Tpetra_INST_COMPLEX_FLOAT AND Teuchos_ENABLE_COMPLEX)

    TRIBITS_ADD_TEST(
      GaleriXpetraDriver
      NAME "Galeri_xpetra_complex_float_Jacobi"
      ARGS "--scalarType=\"complex<float>\" --xml=stratimikos_jacobi.xml"
      NUM_MPI_PROCS 4
      )

  ENDIF()

ENDIF ()


IF (${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND Tpetra_INST_DOUBLE AND Tpetra_INST_INT_INT)
 TRIBITS_ADD_EXECUTABLE(
   tpetra_ifpack2_issue_535
   SOURCES
     test_thyra_tpetra_ifpack2_issue_535.cpp
   COMM serial mpi
   )

 TRIBITS_ADD_TEST(
   tpetra_ifpack2_issue_535
   NAME "tpetra_ifpack2_issue_535"
   NUM_MPI_PROCS 1
   )

ENDIF ()
