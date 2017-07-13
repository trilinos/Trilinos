
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE(
    TogglePFactory_driver
    SOURCES main.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
      ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos)

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_tent_tent_Epetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_tent_tent.xml --linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_sa_tent_Epetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_sa_tent.xml --linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_semi_tent_Epetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_tent.xml --linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_tent_tent_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_tent_tent.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_sa_tent_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_sa_tent.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_semi_tent_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_tent.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_semi_tent_line_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_tent_line.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_semi_sa_line_easy_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_sa_easy.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

    TRIBITS_ADD_TEST(
      TogglePFactory_driver
      NAME "Driver_TogglePFactory_semi_sa_line_easy_Tpetra"
      ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_sa_easy.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )


  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ToggleDriver_cp
    SOURCE_FILES parameters_sa_tent.xml parameters_semi_tent.xml parameters_semi_tent_line.xml parameters_tent_tent.xml parameters_semi_sa_easy.xml
   )
ENDIF()
