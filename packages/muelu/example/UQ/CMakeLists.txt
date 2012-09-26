INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Stokhos AND ${PACKAGE_NAME}_ENABLE_PCE_SUPPORT)
  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    muestok
    SOURCES Simple.cpp
    DEPLIBS muelu xpetra stokhos
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()
