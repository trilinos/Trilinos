
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

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TogglePFactory_test
    SOURCES main.cpp
    ARGS "--matrixType=Laplace3D --nx=40 --ny=40 --nz=60 --mz=1"
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Driver_cp
    SOURCE_FILES parameters_sa_tent.xml parameters_semi_tent.xml paramters_tent_tent.xml
   )
ENDIF()
