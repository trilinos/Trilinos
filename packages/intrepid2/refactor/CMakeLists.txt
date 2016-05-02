ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(unit-test)

## performance tests 
IF (Intrepid2_ENABLE_PerformanceTests)
TRIBITS_ADD_TEST_DIRECTORIES(perf-test)
ENDIF()


## we do not have examples yet
##TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
