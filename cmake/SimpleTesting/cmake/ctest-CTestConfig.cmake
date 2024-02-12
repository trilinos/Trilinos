# This file will be copied to "CTestConfig.cmake" and should be placed
# in the ${CTEST_BINARY_DIRECTORY} directory.
#
# https://cmake.org/cmake/help/latest/module/CTest.html
message(">>> CTestConfig.cmake STARTED")

# Must match what is in CDash project 'Trilinos'
set(CTEST_NIGHTLY_START_TIME "04:00:00 UTC") # 10 PM MDT or 9 PM MST


if (NOT DEFINED CTEST_DROP_METHOD)
    set(CTEST_DROP_METHOD "http")
endif()


# Temporarily changing TEST_DROP_SITE from testing-vm.sandia.gov to
# testing.sandia.gov due to the broken CDash installation
if (CTEST_DROP_METHOD STREQUAL "http")
    set(CTEST_PROJECT_NAME "Trilinos")
    if (NOT DEFINED CTEST_DROP_SITE)
        set(CTEST_DROP_SITE "testing.sandia.gov")
    endif()

    if (${CTEST_DROP_SITE} STREQUAL "testing.sandia.gov")
    	set(CTEST_DROP_LOCATION "/cdash/submit.php?project=${CTEST_PROJECT_NAME}")
    else()
	set(CTEST_DROP_METHOD "https")
	set(CTEST_DROP_LOCATION "/submit.php?project=${CTEST_PROJECT_NAME}")
    endif()

    set(CTEST_TRIGGER_SITE "")
    set(CTEST_DROP_SITE_CDASH TRUE)
    message(">>> CTEST_DROP_SITE     : ${CTEST_DROP_SITE}")
    message(">>> CTEST_DROP_LOCATION : ${CTEST_DROP_LOCATION}")
    message(">>> CTEST_PROJECT_NAME  : ${CTEST_PROJECT_NAME}")
endif()

message(">>> CTestConfig.cmake FINISHED")
