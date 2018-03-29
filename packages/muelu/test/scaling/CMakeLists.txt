# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )

#
# Executable
#


#TODO: how to remove that?
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Galeri)
IF ((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

  TRIBITS_ADD_EXECUTABLE(
    Driver
    SOURCES Driver.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Driver_cp
    SOURCE_FILES scaling.xml scaling-complex.xml scaling-withglobalconstants.xml scaling-complex-withglobalconstants.xml circ_nsp_dependency.xml isorropia.xml iso_poisson.xml conchas_milestone_zoltan.xml conchas_milestone_zoltan2.xml conchas_milestone_zoltan2_complex.xml sa_with_ilu.xml rap.xml smoother.xml smoother_complex.xml tripleMatrixProduct.xml scaling-ml.xml
   )

 TRIBITS_ADD_EXECUTABLE(
   ImportPerformance
   SOURCES ImportPerformance
   COMM mpi
   )

 TRIBITS_ADD_EXECUTABLE(
   TAFCPerformance
   SOURCES TAFCPerformance
   COMM mpi
   )

TRIBITS_ADD_EXECUTABLE(
  ReadMatrix
  SOURCES ReadMatrix.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  MatrixMatrixMultiply
  SOURCES MatrixMatrixMultiply.cpp
  COMM mpi
  )

# not very elegant. The best would probably to make Reuse.cpp work for both Epetra and Tpetra
IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  TRIBITS_ADD_EXECUTABLE(
    Reuse
    SOURCES Reuse.cpp
    COMM serial mpi
  )
  TRIBITS_COPY_FILES_TO_BINARY_DIR(Reuse_cp
    SOURCE_FILES reuse_emin_emin.xml reuse_sa_RP.xml
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(ReadMatrix_cp
  SOURCE_FILES A.mm B.mm simple.xml
)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ImportPerformance_cp
    SOURCE_FILES import.xml
    )

ENDIF()

ENDIF() # Galeri

#
# Tests
#

IF (${PACKAGE_NAME}_ENABLE_Galeri)

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos AND (NOT Xpetra_INT_LONG_LONG))

  IF (${PACKAGE_NAME}_ENABLE_Zoltan)
    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_Defaults"
      ARGS "--linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_isotropic_poisson"
      ARGS "--linAlgebra=Epetra --xml=iso_poisson.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_Milestone"
      ARGS "--linAlgebra=Epetra --xml=conchas_milestone_zoltan.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )

    IF (${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_Isorropia)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverEpetra_IsorropiaPoisson"
        ARGS "--linAlgebra=Epetra --xml=isorropia.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ENDIF()

    IF (${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_Isorropia)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverEpetra_CircNspDependency"
        ARGS "--linAlgebra=Epetra --xml=circ_nsp_dependency.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ENDIF()

  ENDIF()

  TRIBITS_ADD_TEST(
    ImportPerformance
    NAME "ImportPerformance_Epetra"
    ARGS "--linAlgebra=Epetra --nx=30 --ny=30"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "RAPScalingTestEpetra"
    ARGS "--linAlgebra=Epetra --xml=rap.xml --solver=none --nx=50 --ny=50 --rebuild=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Driver
    NAME "SmootherScalingTestEpetra"
    ARGS "--linAlgebra=Epetra --xml=smoother.xml --nx=10 --ny=10 --solver=none"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    ReadMatrix
    NAME "ReadMatrixEpetra"
    ARGS "--linAlgebra=Epetra --xml=simple.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  IF (${PACKAGE_NAME}_ENABLE_Zoltan)
    IF (${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra"
        ARGS "--linAlgebra=Tpetra --xml=scaling-complex.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_WithGlobalConstants"
        ARGS "--linAlgebra=Tpetra --xml=scaling-complex-withglobalconstants.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )
    ELSE()
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra"
        ARGS "--linAlgebra=Tpetra --xml=scaling.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_WithGlobalConstants"
        ARGS "--linAlgebra=Tpetra --xml=scaling-withglobalconstants.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )
    ENDIF()

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetraTripleMatrixProduct"
      ARGS "--linAlgebra=Tpetra --xml=tripleMatrixProduct.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetraILU"
      ARGS "--linAlgebra=Tpetra --xml=sa_with_ilu.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )

  ENDIF()



  IF (${PACKAGE_NAME}_ENABLE_Zoltan2)
    IF (${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_Milestone"
        ARGS "--linAlgebra=Tpetra --xml=conchas_milestone_zoltan2_complex.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )
    ELSE()
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_Milestone"
        ARGS "--linAlgebra=Tpetra --xml=conchas_milestone_zoltan2.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        )
    ENDIF()
  ENDIF()

  TRIBITS_ADD_TEST(
    Driver
    NAME "RAPScalingTestTpetra"
    ARGS "--linAlgebra=Tpetra --xml=rap.xml --solver=none --nx=50 --ny=50 --rebuild=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    TRIBITS_ADD_TEST(
      Driver
      NAME "SmootherScalingTestTpetra"
      ARGS "--linAlgebra=Tpetra --xml=smoother_complex.xml --nx=10 --ny=10 --solver=none"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ELSE()
    TRIBITS_ADD_TEST(
      Driver
      NAME "SmootherScalingTestTpetra"
      ARGS "--linAlgebra=Tpetra --xml=smoother.xml --nx=10 --ny=10 --solver=none"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ENDIF()

  # Reading a matrix saved in real format as complex does not work.
  IF (NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    TRIBITS_ADD_TEST(
      ReadMatrix
      NAME "ReadMatrixTpetra"
      ARGS "--linAlgebra=Tpetra --xml=simple.xml"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ENDIF()

ENDIF()


# Driver to evaluate low level Matrix Matrix multiply kernels potentially using different TPLs
# This differs from the existing Matrix Matrix driver that forms parallel random matrices for testing
# This particular driver does not require MPI and is focused on evaluations using smaller
# MatrixMarket matrices.
IF (${PACKAGE_NAME}_ENABLE_Experimental AND ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor AND ${PACKAGE_NAME}_ENABLE_Tpetra)

  # Reading a matrix saved in real format as complex does not work.
  IF (NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    TRIBITS_ADD_EXECUTABLE(
      MMKernelDriver
      SOURCES MMKernelDriver.cpp
      COMM mpi
      )
  

    TRIBITS_ADD_TEST(
      MMKernelDriver
      NAME "MatrixMatrixKernelDriver"
      ARGS "--linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ENDIF()
ENDIF()


ENDIF() # Galeri
