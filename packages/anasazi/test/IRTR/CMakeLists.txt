
ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ../MVOPTester)

IF(Teuchos_ENABLE_COMPLEX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    IRTR_complex_lap_test
    SOURCES cxx_main_complex_lap.cpp
    ARGS 
      "--verbose --hefty --fakeM --sort=SR"
      "--verbose --skinny --noM --sort=SR"
      "--verbose --hefty --fakeM --sort=SR"
      "--verbose --skinny --noM --sort=SR"
      "--verbose --hefty --fakeM --sort=LR"
      "--verbose --skinny --noM --sort=LR" 
      "--verbose --hefty --fakeM --sort=LR"
      "--verbose --skinny --noM --sort=LR"
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )


ENDIF()
