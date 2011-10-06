INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  )

ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos2
  )

PACKAGE_ADD_EXECUTABLE(
   simple1D
   SOURCES main.cpp
   COMM serial mpi
   DEPLIBS muelu muelu-xpetra muelu-adapters
)

IF (${PACKAGE_NAME}_ENABLE_Epetra)
 IF (${PACKAGE_NAME}_ENABLE_Amesos)

  PACKAGE_ADD_TEST(
    simple1D
    NAME "simple1DEpetra"
    ARGS "--linAlgebra=0"
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  PACKAGE_ADD_TEST(
    simple1D
    NAME "simple1DEpetra"
    ARGS "--linAlgebra=0"
    NUM_MPI_PROCS 4
    COMM mpi
    )

 ENDIF()
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Tpetra)
   IF (${PACKAGE_NAME}_ENABLE_Amesos2)
   
  PACKAGE_ADD_TEST(
    simple1D
    NAME "simple1DTpetra"
    ARGS "--linAlgebra=1"
    NUM_MPI_PROCS 1
    COMM serial mpi
    )

  PACKAGE_ADD_TEST(
    simple1D
    NAME "simple1DTpetra"
    ARGS "--linAlgebra=1"
    NUM_MPI_PROCS 4
    COMM mpi
    )

   ENDIF()
ENDIF()

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     amesos2
#     SOURCES amesos2.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu muelu-xpetra
#   )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     simple1D-belosEpetra1
#     SOURCES main-belosEpetra1.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu muelu-xpetra
#   )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     simple1D-belosXpetraTpetra
#     SOURCES main-belosXpetraTpetra.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu muelu-xpetra
#   )


# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     simple1D-mergedSmoo
#     SOURCES mainMergedSmoo.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu muelu-xpetra
#   )

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#     simple1D-belosXpetraEpetra
#     SOURCES main-belosXpetraEpetra.cpp
#     ARGS ""
#     COMM serial mpi
#     DEPLIBS muelu muelu-xpetra
#   )
