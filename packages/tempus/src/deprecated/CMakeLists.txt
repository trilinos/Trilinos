set(HEADERS "")
set(SOURCES "")

set_and_inc_dirs(DIR ${CMAKE_CURRENT_SOURCE_DIR})
append_glob(HEADERS ${DIR}/*.hpp)
append_glob(SOURCES ${DIR}/*.cpp)

if (NOT ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/deprecated)
  append_set(HEADERS
    Tempus_IntegratorBasicOld_decl.hpp
    Tempus_IntegratorBasicOld_impl.hpp
    )
  append_set(SOURCES
    Tempus_IntegratorBasicOld.cpp
    )
endif()

tribits_add_library(
  tempus
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
