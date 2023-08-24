include_guard()
message("+--------------------------------------+")
message("| ctest-cdash-setup.cmake START        |")
message("+--------------------------------------+")
include(${CMAKE_CURRENT_LIST_DIR}/ctest-functions.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ctest-common.cmake)



# ===============================================================================
# --
# -- Set up the Dashboard and Testing Infrastructure
# --
# ===============================================================================
print_options_list()

# -----------------------------------------------------------
# -- Specify the Generator
# -----------------------------------------------------------
set(CTEST_CMAKE_GENERATOR "Ninja")


set(CTEST_CONFIGURE_COMMAND_ARGS
    "${CMAKE_COMMAND}"
    "-C \"${configure_file}\""
    "-C \"${package_enables_file}\""
    "-G \"${CTEST_CMAKE_GENERATOR}\""
    "${CTEST_SOURCE_DIRECTORY}"
)

list(JOIN CTEST_CONFIGURE_COMMAND_ARGS " " CTEST_CONFIGURE_COMMAND)

banner("CMake Configuration Command")
message(">>> CTEST_CONFIGURE_COMMAND : ${CTEST_CONFIGURE_COMMAND}")
message(">>>")



# -----------------------------------------------------------
# -- Configure Preparation
# -----------------------------------------------------------

# Optionally skip cleaning the build directory
if(NOT skip_clean_build_dir)
    ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

# Optionally upload the config files
# TODO: Note how this works / what it's doing in CMake-land.
message(">>> Write `configure_command_file`:")
if(skip_upload_config_files)
    message(">>> - SKIPPED")
else()
    message(">>> - WRITTEN")
    #message(">>> - configure_command_file : ${configure_command_file}")
    #message(">>> - CTEST_CONFIGURE_COMMAND: ${CTEST_CONFIGURE_COMMAND}")
    file(WRITE ${configure_command_file} ${CTEST_CONFIGURE_COMMAND})
endif()
message(">>>")


# -----------------------------------------------------------
# -- COPY files
# -- ctest-CTestConfig.cmake -> ${CTEST_BINARY_DIRECTORY}/CTestConfig.cmake
# -- See: https://cmake.org/cmake/help/latest/command/configure_file.html
# -----------------------------------------------------------
banner("Copy CTestConfig.cmake to BUILD dir")
configure_file(${CMAKE_CURRENT_LIST_DIR}/ctest-CTestConfig.cmake ${CTEST_BINARY_DIRECTORY}/CTestConfig.cmake COPYONLY)


# -----------------------------------------------------------
# -- Start the testing engine
# -- CTEST_DROP_LOCATION defined in ctest-CTestConfig.cmake
# -- See: https://cmake.org/cmake/help/latest/command/ctest_start.html
# -----------------------------------------------------------
banner("ctest_start() START")
ctest_start(${dashboard_model} GROUP ${dashboard_track})
banner("ctest_start() FINISH")
message(">>> CTEST_DROP_LOCATION = ${CTEST_DROP_LOCATION}")
message(">>>")

# -----------------------------------------------------------
# -- Configure CDash Drop URLs
# -- CTEST_DROP_SITE is defined in ctest-CTestConfig.cmake
# -- CTEST_PROJECT_NAME is defined in ctest-CTestConfig.cmake
# -- CTEST_BUILD_NAME is defined in ctest-common.cmake
# -- URL_location is defined here
# -- build_stamp is defined here
# -- machine_name is defined here
# -----------------------------------------------------------

string(FIND ${CTEST_DROP_LOCATION} submit.php index)
string(SUBSTRING ${CTEST_DROP_LOCATION} 0 ${index} URL_location)

file(STRINGS ${CTEST_BINARY_DIRECTORY}/Testing/TAG tag_strings LIMIT_COUNT 2)
foreach(item ${tag_strings})
    set(build_stamp_tmp "${build_stamp_tmp}-${item}")
endforeach()
# strip initial "-" using this method to avoid having to calculate string length
string(SUBSTRING ${build_stamp_tmp} 1 1024 build_stamp)

generate_build_url1(build_url1 ${CTEST_DROP_SITE} ${URL_location} ${CTEST_PROJECT_NAME} ${CTEST_BUILD_NAME} ${build_stamp} ${machine_name})
generate_build_url2(build_url2 ${CTEST_DROP_SITE} ${URL_location} ${CTEST_PROJECT_NAME} ${CTEST_BUILD_NAME} ${build_stamp})
generate_build_url3(build_url3 ${CTEST_DROP_SITE} ${URL_location} ${CTEST_PROJECT_NAME} ${CTEST_BUILD_NAME} ${build_stamp})
message(">>> CDash URL1 = ${build_url1}")
message(">>> CDash URL2 = ${build_url2}")
message(">>> CDash URL3 = ${build_url3}")

# -----------------------------------------------------------
# -- Optionally update the repository
# -- skip_update_step defined in ctest-common.cmake
# -- skip_by_parts_submit defined in ctest-common.cmake
# -----------------------------------------------------------
# TODO: Do we really do this in the Trilinos CTest framework?
#       Or is this all handled by the framework code around
#       pulling in the source/target branches and merging them?
#       If we don't use this, consider removing.
banner("Update repository START")
if(NOT skip_update_step)
    message(">>> Updating repository")
    ctest_update(RETURN_VALUE update_error)
    if(${update_error} EQUAL -1)
        message(WARNING ">>> Update failed. ")
    else()
        message(">>> Updated ${update_error} files.")
    endif()
    submit_by_parts("Update")
else()
    set(CTEST_UPDATE_VERSION_ONLY ON)
    ctest_update(RETURN_VALUE update_error)
    if(${update_error} EQUAL -1)
        message(WARNING ">>> Update failed. ")
    else()
        message(">>> Updated ${update_error} files.")
    endif()
    submit_by_parts("Update")
endif()
banner("Update repository FINISH")

message("+--------------------------------------+")
message("| ctest-cdash-setup.cmake FINISH       |")
message("+--------------------------------------+")
