INCLUDE(PackageAddExecutableAndTest)

SET(SOURCES "main.cpp")
FILE(GLOB SOURCES "*.cpp")

PACKAGE_ADD_EXECUTABLE(
  XpetraCompilationTests
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS muelu-xpetra
  )
