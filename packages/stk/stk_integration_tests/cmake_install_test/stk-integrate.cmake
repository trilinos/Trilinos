message(STATUS "Integration build of STK in Trilinos")

if(NOT DEFINED CONFIGURE_ARGS)
  set(CONFIGURE_ARGS "")
endif()

set($ENV{LC_MESSAGES} "en_EN")
message(STATUS ${DASHBOARD})

if(NOT DEFINED(CUDA))
  set(CUDA OFF)
endif()
if(NOT DEFINED(CLEAR_CACHE))
  set(CLEAR_CACHE ON)
endif()
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE release)
endif()
if(NOT DEFINED NUM_JOBS)
  set(NUM_JOBS 12)
endif()
if(NOT DEFINED CDASH_TRACK)
  set(CDASH_TRACK "foobar")
endif()

set(STK_SCRIPT_DIR "${SIERRA_PROJ}/stk/stk_integration_tests/cmake_install_test")
set(TRILINOS_DIR "${OUTPUT_DIR}/Trilinos")

set(CTEST_BINARY_DIRECTORY "${OUTPUT_DIR}/build")
set(CTEST_SOURCE_DIRECTORY "${TRILINOS_DIR}")


set(CTEST_UPDATE_COMMAND "${STK_SCRIPT_DIR}/create_workspace.sh")
set(CTEST_UPDATE_OPTIONS "${TRILINOS_DIR} ${SIERRA_PROJ}")
set(CTEST_CHECKOUT_COMMAND "${CTEST_UPDATE_COMMAND} ${CTEST_UPDATE_OPTIONS}")


set(CTEST_CONFIGURE_COMMAND "${STK_SCRIPT_DIR}/run_cmake_stk")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_TEST_TIMEOUT "600")

message(STATUS " -- Start dashboard --")
ctest_start("Experimental" GROUP ${CDASH_TRACK})
ctest_configure()
ctest_build(PARALLEL_LEVEL ${NUM_JOBS} RETURN_VALUE ESTAT)
ctest_test(TEST_LOAD ${NUM_JOBS})

message(STATUS " -- Finished --")
if(NOT ${ESTAT} EQUAL 0)
  message(FATAL_ERROR "Testing failed")
endif()
