if ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
  message(FATAL_ERROR "Error, CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    " == CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
    " TrilinosDriver does not support in source builds!" )
endif()

cmake_minimum_required(VERSION 2.7.20090924)
project(TrilinosDriver NONE)

# We depend on some of the cmake modules in the parent Trilinos source tree:
#
get_filename_component(TRILINOS_HOME_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
get_filename_component(TD_BASE_DIR
  "${TRILINOS_HOME_DIR}" PATH)

set(CMAKE_MODULE_PATH
  ${TRILINOS_HOME_DIR}
  ${TRILINOS_HOME_DIR}/cmake
  ${TRILINOS_HOME_DIR}/cmake/utils
  ${TRILINOS_HOME_DIR}/cmake/package_arch
  ${TRILINOS_HOME_DIR}/cmake/config_tests
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

include(CTest)
include(TrilinosDriverSupport)

# Now make exactly one add_subdirectory call based on the site name of the
# machine we're running on...
#
# By default, use the exact site name as the subdir value:
#
site_name(site)
set(subdir "${site}")
message("site='${site}'")

# But if that directory does not exist as named, and there's a regex match
# with the name of a subdirectory, use the exact subdirectory name instead:
#
file(GLOB filesAndDirs RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*")

set(dirs "")
foreach(dir ${filesAndDirs})
  if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
    set(dirs ${dirs} "${dir}")
  endif()
endforeach()

foreach(dir ${dirs})
  if("${site}" MATCHES "${dir}")
    set(subdir "${dir}")
    message("site='${site}' MATCHES directory name dir='${dir}'")
  endif()
endforeach()

# The one add_subdirectory call:
#
message("subdir='${subdir}'")
add_subdirectory("${subdir}")
