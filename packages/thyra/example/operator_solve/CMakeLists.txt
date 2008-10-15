INCLUDE(Trilinos_Add_Executable_And_Test)


TRILINOS_ADD_EXECUTABLE_AND_TEST(
  silly1DFFT_serial
  SOURCES silly1DFFT_serial.cpp
  ARGS --verbose
  COMM serial
  )
