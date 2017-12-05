ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

IF(Stokhos_ENABLE_Sacado)
  IF(Stokhos_ENABLE_Ensemble_Scalar_Type)
    IF (Stokhos_ENABLE_TpetraKernels)
      IF (Stokhos_ENABLE_Tpetra)
        IF (HAVE_TPETRA_SERIAL)
          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            ThyraMP_Test
            SOURCES main.cpp
            TESTONLYLIBS thyra_nonlinear_test_models
            COMM serial mpi
            NUM_MPI_PROCS 4
            )
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()
ENDIF()

