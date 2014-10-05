INCLUDE(TribitsAddExecutableAndTest)

SET(SOURCES "main.cpp")
FILE(GLOB SOURCES "*.cpp")

TRIBITS_ADD_EXECUTABLE(
  XpetraCompilationTests
  SOURCES ${SOURCES}
  COMM serial mpi
  DEPLIBS xpetra
  )
