
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/../thyra/example/operator_solve)


TRILINOS_ADD_EXECUTABLE(
  test_single_stratimikos_solver_driver
  SOURCES
    test_single_stratimikos_solver_driver.cpp
    test_single_stratimikos_solver.hpp
    test_single_stratimikos_solver.cpp
  COMM serial mpi
  )


COPY_FILES_TO_BINARY_DIR(CopyTestFourByFourMtx
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS test_single_stratimikos_solver_driver
  )


COPY_FILES_TO_BINARY_DIR(CopyTestXMLFiles
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
  SOURCE_PREFIX "_"
  EXEDEPS test_single_stratimikos_solver_driver
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Amesos)
  TRILINOS_ADD_TEST(
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


  TRILINOS_ADD_TEST(
    test_single_stratimikos_solver_driver
    NAME test_single_stratimikos_solver_driver_aztecoo_np
    ARGS "--input-file=FourByFour.aztecoo.np.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    TRILINOS_ADD_TEST(
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

    TRILINOS_ADD_TEST(
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


  TRILINOS_ADD_TEST(
    test_single_stratimikos_solver_driver
    NAME test_single_stratimikos_solver_driver_belos_np
    ARGS "--input-file=FourByFour.belos.np.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    TRILINOS_ADD_TEST(
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

    TRILINOS_ADD_TEST(
      test_single_stratimikos_solver_driver
      NAME test_single_stratimikos_solver_driver_belos_ml
      ARGS
        "--input-file=FourByFour.belos.ml.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

  ENDIF()


ENDIF()
