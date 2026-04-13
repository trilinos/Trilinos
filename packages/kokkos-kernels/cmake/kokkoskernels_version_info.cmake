# https://jonathanhamberg.com/post/cmake-embedding-git-hash/

find_package(Git QUIET)

set(CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR})
set(pre_configure_file ${CURRENT_LIST_DIR}/KokkosKernels_Version_Info.hpp.in)
set(post_configure_file ${CMAKE_BINARY_DIR}/KokkosKernels_Version_Info.hpp)

function(check_git_write git_hash git_clean_status)
  file(WRITE ${CMAKE_BINARY_DIR}/git-state.txt "${git_hash}-${git_clean_status}")
endfunction()

function(check_git_read git_hash)
  if(EXISTS ${CMAKE_BINARY_DIR}/git-state.txt)
    file(STRINGS ${CMAKE_BINARY_DIR}/git-state.txt CONTENT)
    list(GET CONTENT 0 var)

    message(DEBUG "Cached Git hash: ${var}")
    set(${git_hash} ${var} PARENT_SCOPE)
  else()
    set(${git_hash} "INVALID" PARENT_SCOPE)
  endif()
endfunction()

function(check_git_version)
  if(NOT Git_FOUND OR NOT EXISTS ${KOKKOSKERNELS_TOP_SOURCE_DIR}/.git)
    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    return()
  endif()

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

  if(IS_DIRTY EQUAL 0)
    set(GIT_CLEAN_STATUS "CLEAN")
  else()
    set(GIT_CLEAN_STATUS "DIRTY")
  endif()

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
    WORKING_DIRECTORY ${KOKKOSKERNELS_TOP_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  check_git_read(GIT_HASH_CACHE)

  # Only update the version header if the hash has changed. This will
  # prevent us from rebuilding the project more than we need to.
  if(NOT "${GIT_COMMIT_HASH}-${GIT_CLEAN_STATUS}" STREQUAL ${GIT_HASH_CACHE} OR NOT EXISTS ${post_configure_file})
    # Set the GIT_HASH_CACHE variable so the next build won't have
    # to regenerate the source file.
    check_git_write(${GIT_COMMIT_HASH} ${GIT_CLEAN_STATUS})

    configure_file(${pre_configure_file} ${post_configure_file} @ONLY)
    message(STATUS "Configured git information in ${post_configure_file}")
  endif()
endfunction()

# Pass BENCHMARK_VERSION variable to configure benchmark library version
function(check_version_info)
  add_custom_target(AlwaysCheckGit
    COMMAND
      ${CMAKE_COMMAND} -DRUN_CHECK_GIT_VERSION=1
      -DKOKKOSKERNELS_TOP_SOURCE_DIR=${KOKKOSKERNELS_TOP_SOURCE_DIR}
      -DBENCHMARK_VERSION=${BENCHMARK_VERSION} -P
      ${CURRENT_LIST_DIR}/kokkoskernels_version_info.cmake
    BYPRODUCTS ${post_configure_file})

  add_dependencies(kokkoskernels AlwaysCheckGit)
  check_git_version()
endfunction()

# This is used to run this function from an external cmake process.
if(RUN_CHECK_GIT_VERSION)
  check_git_version()
endif()
