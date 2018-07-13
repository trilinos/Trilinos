TRIBITS_SUBPACKAGE(Exodus)

if (${CMAKE_PROJECT_NAME} STREQUAL "SEACASProj")
INCLUDE(InstallSymLink)
endif()

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

FILE(GLOB SOURCES src/ex_*.c)

if (${SEACAS_HIDE_DEPRECATED_CODE})
else()
  FILE(GLOB DEP_SOURCES src/deprecated/ex_*.c)
endif()
TRIBITS_CONFIGURE_FILE(exodus_config.h)

if (NOT ${EXODUS_THREADSAFE})
   list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/ex_threadsafe.c)
endif()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  include/exodusII.h
  include/exodusII_par.h
  ${CMAKE_CURRENT_BINARY_DIR}/exodus_config.h
)

TRIBITS_ADD_LIBRARY(
  exodus
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${DEP_SOURCES}
)

if (${TPL_Netcdf_PARALLEL})
   add_definitions(-DPARALLEL_AWARE_EXODUS)
endif()

if (${CMAKE_PROJECT_NAME} STREQUAL "SEACASProj")
   # Many users of the SEACAS installations still want a static exodus library...
   if (BUILD_SHARED_LIBS)
      ADD_LIBRARY(
              exodus_static
              STATIC
              ${SOURCES} ${DEP_SOURCES}
      )
      set_target_properties(exodus_static PROPERTIES OUTPUT_NAME exodus)
      INSTALL(TARGETS exodus_static DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   endif(BUILD_SHARED_LIBS)

   # For backwards compatability, the C exodus library was called libexoIIv2c.a
   # Create a symbolic link from libexodus.a to libexoIIv2c.a -- do only for static lib
   InstallSymLink(libexodus.a ${CMAKE_INSTALL_PREFIX}/lib/libexoIIv2c.a)
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
