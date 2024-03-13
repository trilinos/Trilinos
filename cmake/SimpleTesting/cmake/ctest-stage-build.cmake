message("+--------------------------------------+")
message("| ctest-stage-build.cmake START        |")
message("+--------------------------------------+")
set(STAGE_BUILD_ERROR OFF)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-cdash-setup.cmake)


# -----------------------------------------------------------
# -- Build
# -----------------------------------------------------------
banner("START build step")

# Execute the build step
ctest_build(RETURN_VALUE build_error)

# Print out final stage banner
if(${build_error} EQUAL 0)
    banner("END build step - SUCCESS")
else()
    message(WARNING "Build failed with error ${build_error}")
    banner("END build step - FAILURE")
    set(STAGE_BUILD_ERROR ON)
endif()

# Submit to CDash
submit_by_parts( "Build" )

message("+--------------------------------------+")
message("| ctest-stage-build.cmake FINISH       |")
message("+--------------------------------------+")
