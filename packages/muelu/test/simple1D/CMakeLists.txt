INCLUDE(PackageAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    simple1D
    SOURCES main.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu muelu-cthulhu
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    simple1D-belosEpetra1
    SOURCES main-belosEpetra1.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu muelu-cthulhu
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
    simple1D-belosCthulhuTpetra
    SOURCES main-belosCthulhuTpetra.cpp
    ARGS ""
    COMM serial mpi
    DEPLIBS muelu muelu-cthulhu
  )

