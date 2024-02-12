# https://jonathanhamberg.com/post/cmake-embedding-git-hash/

find_package(Git QUIET)

SET(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
SET(pre_configure_file ${CURRENT_LIST_DIR}/KokkosKernels_Version_Info.hpp.in)
SET(post_configure_file ${CMAKE_BINARY_DIR}/KokkosKernels_Version_Info.hpp)

FUNCTION(check_git_write git_hash git_clean_status)
  FILE(
    WRITE
    ${CMAKE_BINARY_DIR}/git-state.txt
    "${git_hash}-${git_clean_status}")
ENDFUNCTION()

FUNCTION(check_git_read git_hash)
  IF(EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
    FILE(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
    LIST(GET CONTENT 0 var)

    message(DEBUG "Cached Git hash: ${var}")
    SET(${git_hash} ${var} PARENT_SCOPE)
  else()
    SET(${git_hash} "INVALID" PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

FUNCTION(check_git_version)
  IF(NOT Git_FOUND OR NOT EXISTS ${KOKKOSKERNELS_TOP_SOURCE_DIR}/.git)
    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    return()
  ENDIF()

  # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest commit description
  execute_process(
    COMMAND ${GIT_EXECUTABLE} show -s --format=%s
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DESCRIPTION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest commit date
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%cI
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Check if repo is dirty / clean
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    RESULT_VARIABLE IS_DIRTY
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  IF(IS_DIRTY EQUAL 0)
    SET(GIT_CLEAN_STATUS "CLEAN")
  else()
    SET(GIT_CLEAN_STATUS "DIRTY")
  ENDIF()

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  check_git_read(GIT_HASH_CACHE)

  # Only update the version header if the hash has changed. This will
  # prevent us from rebuilding the project more than we need to.
  IF(NOT "${GIT_COMMIT_HASH}-${GIT_CLEAN_STATUS}" STREQUAL ${GIT_HASH_CACHE}
    OR NOT EXISTS ${post_configure_file})
    # Set the GIT_HASH_CACHE variable so the next build won't have
    # to regenerate the source file.
    check_git_write(${GIT_COMMIT_HASH} ${GIT_CLEAN_STATUS})

    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    message(STATUS "Configured git information in ${post_configure_file}")
  ENDIF()
ENDFUNCTION()

# Pass BENCHMARK_VERSION variable to configure benchmark library version
FUNCTION(check_version_info)
  add_custom_target(
    AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
    -DRUN_CHECK_GIT_VERSION=1
    -DKOKKOSKERNELS_TOP_SOURCE_DIR=${KOKKOSKERNELS_TOP_SOURCE_DIR}
    -DBENCHMARK_VERSION=${BENCHMARK_VERSION}
    -P ${CURRENT_LIST_DIR}/kokkoskernels_version_info.cmake
    BYPRODUCTS ${post_configure_file})

  add_dependencies(kokkoskernels AlwaysCheckGit)
  check_git_version()
ENDFUNCTION()

# This is used to run this function from an external cmake process.
IF(RUN_CHECK_GIT_VERSION)
  check_git_version()
ENDIF()
