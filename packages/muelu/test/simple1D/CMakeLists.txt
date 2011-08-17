INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../utils)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    simple1D
    SOURCES main.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu muelu-xpetra
  )

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
