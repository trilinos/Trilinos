
INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  aztec_simple
  SOURCES ml_aztec_simple.c place_holder.cpp
  ARGS 
  COMM serial
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  operator
  SOURCES ml_operator.cpp
  ARGS 
  COMM serial
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  preconditioner
  SOURCES ml_preconditioner.cpp
  ARGS 
  COMM serial
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  readex
  SOURCES ml_readex.c place_holder.cpp
  ARGS 
  COMM serial
  )

TRILINOS_ADD_EXECUTABLE_AND_TEST(
  user_smoothing
  SOURCES ml_user_smoothing.cpp
  ARGS 
  COMM serial
  )

