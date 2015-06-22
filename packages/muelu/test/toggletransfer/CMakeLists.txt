
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

  TRIBITS_ADD_TEST(
    TogglePFactory_driver
    NAME "Driver_TogglePFactory_tent_tent"
    ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_tent_tent.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    TogglePFactory_driver
    NAME "Driver_TogglePFactory_sa_tent"
    ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_sa_tent.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    TogglePFactory_driver
    NAME "Driver_TogglePFactory_semi_tent"
    ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1 --xmlFile=parameters_semi_tent.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ToggleDriver_cp
    SOURCE_FILES parameters_sa_tent.xml parameters_semi_tent.xml parameters_tent_tent.xml
   )
ENDIF()
