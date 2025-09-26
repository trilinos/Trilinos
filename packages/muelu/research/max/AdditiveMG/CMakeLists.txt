TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(SOURCES ${DIR}/*.cpp)

#
# Executables
#

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiplicativeStride
    SOURCES MultiplicativeStride.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MultiplicativeBricks
    SOURCES MultiplicativeBricks.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SmoothedAdditiveStride
    SOURCES SmoothedAdditiveStride.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SmoothedAdditiveBricks
    SOURCES SmoothedAdditiveBricks.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    tentative
    SOURCES tentative.cpp
    NUM_MPI_PROCS 4
    COMM serial mpi
  )
  #TRIBITS_COPY_FILES_TO_BINARY_DIR(Max_cp
  #  SOURCE_FILES max.xml
  # )

ENDIF()
