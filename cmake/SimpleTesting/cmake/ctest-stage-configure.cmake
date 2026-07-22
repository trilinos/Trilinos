message("+--------------------------------------+")
message("| ctest-stage-configure.cmake START    |")
message("+--------------------------------------+")
set(STAGE_CONFIGURE_ERROR OFF)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)

# -----------------------------------------------------------
# -- Configure
# -----------------------------------------------------------
banner("START configure step")

message(">>> CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(">>> CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

ctest_configure(SOURCE ${CTEST_SOURCE_DIRECTORY}
                BUILD  ${CTEST_BINARY_DIRECTORY}
                RETURN_VALUE configure_error
                CAPTURE_CMAKE_ERROR tmp_cmake_error
               )

ctest_read_custom_files(BUILD "${CTEST_BINARY_DIRECTORY}")
# NOTE: The CTestCustom.cmake file read in by the above command is
# automatically written in the binary dir by the configure of Trilinos in
# ctest_configure() above.

message(">>> configure_error: ${configure_error}")
message(">>> tmp_cmake_error: ${tmp_cmake_error}")   # NEW (testing this out)

if( (${configure_error} EQUAL 0) AND (${tmp_cmake_error} EQUAL 0) )
    banner("END configure step - SUCCESS")
else()
    banner("END configure step - FAILURE")
    set(STAGE_CONFIGURE_ERROR ON)
endif()

submit_by_parts("Configure")

message("+--------------------------------------+")
message("| ctest-stage-configure.cmake FINISH   |")
message("+--------------------------------------+")
