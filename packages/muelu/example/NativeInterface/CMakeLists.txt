INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    PlainAggregation
    SOURCES PlainAggregation.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    SmoothedAggregation
    SOURCES SmoothedAggregation.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
  TRIBITS_ADD_EXECUTABLE(
    SeparateSetups
    SOURCES SeparateSetups.cpp
    COMM serial mpi
    DEPLIBS muelu xpetra
    )
ENDIF()
