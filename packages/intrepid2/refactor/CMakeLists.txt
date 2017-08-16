ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(unit-test)

## performance tests 
TRIBITS_ADD_EXAMPLE_DIRECTORIES(perf-test)


## we do not have examples yet
##TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
