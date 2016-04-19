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

#
# Tests
#

IF (${PACKAGE_NAME}_ENABLE_Galeri)

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  IF (${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_AztecOO)
    TRIBITS_ADD_EXECUTABLE(
      Viz
      SOURCES Viz.cpp
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      Viz
      NAME "Viz3DTpetra"
      ARGS "--matrixType=Laplace3D --nx=20 --ny=20 --nz=20 --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      Viz
      NAME "Viz2DTpetra"
      ARGS "--matrixType=Laplace2D --nx=40 --ny=40 --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
      
    # Test Epetra variant, too.
    IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack  AND ${PACKAGE_NAME}_ENABLE_Amesos AND (NOT Xpetra_INT_LONG_LONG))

      TRIBITS_ADD_TEST(
        Viz
        NAME "Viz3DEpetra"
        ARGS "--matrixType=Laplace3D --nx=20 --ny=20 --nz=20 --linAlgebra=Epetra"
        NUM_MPI_PROCS 4
        COMM serial mpi
        )
  
      TRIBITS_ADD_TEST(
        Viz
        NAME "Viz2DEpetra"
        ARGS "--matrixType=Laplace2D --nx=40 --ny=40 --linAlgebra=Epetra"
        NUM_MPI_PROCS 4
        COMM serial mpi
        )

    ENDIF() #endif Epetra
      
  ENDIF()

ENDIF()



TRIBITS_COPY_FILES_TO_BINARY_DIR(Viz_cp
    SOURCE_FILES viztest.xml
    )


ENDIF() # Galeri
