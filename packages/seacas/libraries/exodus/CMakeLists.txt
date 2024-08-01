TRIBITS_SUBPACKAGE(Exodus)

if (SEACASExodus_ENABLE_STATIC)
INCLUDE(InstallSymLink)
endif()

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

FILE(GLOB SOURCES src/ex_*.c)

if (NOT ${PACKAGE_NAME}_HIDE_DEPRECATED_CODE)
  FILE(GLOB DEP_SOURCES src/deprecated/ex_*.c)
endif()
TRIBITS_CONFIGURE_FILE(exodus_config.h)

if (NOT ${EXODUS_THREADSAFE})
   list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/ex_threadsafe.c)
endif()

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  include/exodusII.h
  include/exodusII_int.h
  include/exodusII_par.h
  ${CMAKE_CURRENT_BINARY_DIR}/exodus_config.h
)

TRIBITS_ADD_LIBRARY(
  exodus
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${DEP_SOURCES}
)

set_property(TARGET exodus PROPERTY C_STANDARD 99)
set_property(TARGET exodus PROPERTY C_EXTENSIONS ON)

if (SEACASExodus_ENABLE_SHARED)
   if (NOT BUILD_SHARED_LIBS)
      ADD_LIBRARY(
              exodus_shared
              SHARED
              ${SOURCES} ${DEP_SOURCES}
      )
      target_include_directories(exodus_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
      target_link_libraries(exodus_shared PUBLIC Netcdf::all_libs)
      set_property(TARGET exodus_shared PROPERTY C_STANDARD 99)
      # This keeps the library out of the `all_libs` targets...
      set_target_properties(exodus_shared PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
      set_target_properties(exodus_shared PROPERTIES OUTPUT_NAME exodus)
      INSTALL(TARGETS exodus_shared DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   endif()
endif()

# Many users of the SEACAS installations still want a static exodus library...
if (SEACASExodus_ENABLE_STATIC)
   # NOTE: If use `TRIBITS_ADD_LIBRARY` here, then both exodus and exodus_static linked to all executables.
   ADD_LIBRARY(
           exodus_static
           STATIC
           ${SOURCES} ${DEP_SOURCES}
   )
   target_include_directories(exodus_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
   target_link_libraries(exodus_static PUBLIC Netcdf::all_libs)
   set_property(TARGET exodus_static PROPERTY C_STANDARD 99)
   # This keeps the library out of the `all_libs` targets...
   set_target_properties(exodus_static PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
   set_target_properties(exodus_static PROPERTIES OUTPUT_NAME exodus)
   INSTALL(TARGETS exodus_static DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

   # For backwards compatibility, the C exodus library was called libexoIIv2c.a
   # Create a symbolic link from libexodus.a to libexoIIv2c.a -- do only for static lib
   InstallSymLink(libexodus.a ${CMAKE_INSTALL_PREFIX}/lib/libexoIIv2c.a)
endif()

if (${CMAKE_PROJECT_NAME} STREQUAL "Seacas")
   # add a target to generate API documentation with Doxygen
   if(Seacas_ENABLE_DOXYGEN)
     find_package(Doxygen)
     if(DOXYGEN_FOUND)
       add_custom_command(TARGET exodus POST_BUILD
          COMMAND ${DOXYGEN_EXECUTABLE} ARGS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating Exodus API documentation with Doxygen" VERBATIM
          )
     endif(DOXYGEN_FOUND)
   endif()
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
