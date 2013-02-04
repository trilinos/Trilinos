INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    QRdecomposition
    SOURCES QRdecomposition.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
  )
ENDIF()