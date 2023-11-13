set(HEADERS "")
set(SOURCES "")

tribits_set_and_inc_dirs(DIR ${CMAKE_CURRENT_SOURCE_DIR})
append_glob(HEADERS ${DIR}/*.hpp)
append_glob(SOURCES ${DIR}/*.cpp)

if (NOT ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
  tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR}/deprecated)
  append_set(HEADERS
    )
  append_set(SOURCES
    )
endif()

tribits_add_library(
  tempus
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
