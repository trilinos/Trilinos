INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(SOURCES ${DIR}/*.cpp)

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2)

  TRIBITS_ADD_EXECUTABLE(
    hhg_driver 
    SOURCES hhg_driver.cpp
  )

#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    Test_xpetra
#    SOURCES Test_xpetra.cpp
#    NUM_MPI_PROCS 4
#    COMM serial mpi
#    )


ENDIF()
