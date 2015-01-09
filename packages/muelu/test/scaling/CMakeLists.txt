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

# ScalingTest is compile with long long int when Teuchos_ENABLE_LONG_LONG_INT=ON
# The name of the executable is changed to make sure that long int support is correctly
# set when we do scalability tests.
IF(Teuchos_ENABLE_LONG_LONG_INT)
  IF(NOT ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION) #FIXME
    #FIXME no Tpetra instantiation for long long int
    #FIXME we need a HAVE_MUELU_LONG_LONG_INT option
    SET(MueLu_ScalingTestName ScalingTest-LongLongInt)
  ELSE()
    SET(MueLu_ScalingTestName ScalingTest-Int)
  ENDIF()
ELSE()
  SET(MueLu_ScalingTestName ScalingTest-Int)
ENDIF()
#TODO: It is actually possible to build and test ScalingTest-Int also when LONG_LONG_INT=ON

IF (${PACKAGE_NAME}_ENABLE_Galeri)
IF ((${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2) OR
    (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos))

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

  TRIBITS_ADD_EXECUTABLE(
    ${MueLu_ScalingTestName}  #TODO: Suffix Int/LongInt when done in the source code
    SOURCES ScalingTest.cpp
    COMM serial mpi
    )


  TRIBITS_ADD_EXECUTABLE(
    Driver
    SOURCES Driver.cpp
    COMM serial mpi
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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Driver_cp
    SOURCE_FILES scaling.xml isorropia.xml iso_poisson.xml conchas_milestone_zoltan.xml conchas_milestone_zoltan2.xml sa_with_ilu.xml rap.xml smoother.xml
   )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ReadMatrix_cp
    SOURCE_FILES A.mm simple.xml
  )

ENDIF()

ENDIF() # Galeri

#
# Tests
#

INCLUDE(XpetraAddTest)

IF (${PACKAGE_NAME}_ENABLE_Galeri)

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos AND (NOT Teuchos_ENABLE_LONG_LONG_INT))
  XPETRA_ADD_TEST_EPETRA(${MueLu_ScalingTestName} 4)

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

  ENDIF()

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
  XPETRA_ADD_TEST_TPETRA(${MueLu_ScalingTestName} 4)

  IF (${PACKAGE_NAME}_ENABLE_Zoltan)
    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetra"
      ARGS "--linAlgebra=Tpetra"
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
    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetra_Milestone"
      ARGS "--linAlgebra=Tpetra --xml=conchas_milestone_zoltan2.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      )
  ENDIF()

  TRIBITS_ADD_TEST(
    Driver
    NAME "RAPScalingTestTpetra"
    ARGS "--linAlgebra=Tpetra --xml=rap.xml --solver=none --nx=50 --ny=50 --rebuild=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Driver
    NAME "SmootherScalingTestTpetra"
    ARGS "--linAlgebra=Tpetra --xml=smoother.xml --nx=10 --ny=10 --solver=none"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )


  TRIBITS_ADD_TEST(
    ReadMatrix
    NAME "ReadMatrixTpetra"
    ARGS "--linAlgebra=Tpetra --xml=simple.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()

ENDIF() # Galeri
