TRIBITS_ADD_EXECUTABLE_AND_TEST(
  klu2_simple
  HEADERS ${HEADERS}
  SOURCES klu2_simple.cpp
  COMM serial
  )
