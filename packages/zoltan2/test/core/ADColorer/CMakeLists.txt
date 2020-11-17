
TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ADColorer
      NAME ADColorer_simple
      SOURCES ADColorer.cpp
      ARGS "--file=simple.mtx"
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "PASS"
      FAIL_REGULAR_EXPRESSION "FAIL"
    )

IF (${PROJECT_NAME}_ENABLE_Galeri)
TRIBITS_ADD_TEST(
    ADColorer
    NAME ADColorer_galeri
    ARGS "--xdim=10 --ydim=8 --zdim=6"
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAIL"
)
ENDIF()
  
