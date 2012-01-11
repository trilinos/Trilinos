
INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_LIBRARY(
  teuchos_algorithma
  HEADERS AlgorithmA.hpp
  SOURCES AlgorithmA.cpp
  DEPLIBS teuchos
  TESTONLY
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FancyOutputting_test  
  SOURCES FancyOutputting_test.cpp
  DEPLIBS teuchos_algorithma
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  XHOSTTYPE Windows Windows_NT
  )
