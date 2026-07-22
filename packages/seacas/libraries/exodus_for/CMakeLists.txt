INCLUDE(FortranSettings)

if (SEACASExodus_ENABLE_STATIC)
INCLUDE(InstallSymLink)
endif()

TRIBITS_SUBPACKAGE(Exodus_for)

TRIBITS_INCLUDE_DIRECTORIES(
   "${CMAKE_CURRENT_SOURCE_DIR}/include"
   "${CMAKE_CURRENT_SOURCE_DIR}/test"
   )

SET(SOURCES "")
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(SOURCES ${DIR}/src/*.c ${DIR}/src/*.F)
LIST(REMOVE_ITEM SOURCES "${DIR}/src/exo_jack-windows.c")

SET(HEADERS include/exodusII.inc )

TRIBITS_ADD_LIBRARY(
  exodus_for
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# Many users of the SEACAS installations still want a static exodus library...
if (SEACASExodus_ENABLE_STATIC)
   ADD_LIBRARY(
           exodus_for_static
           STATIC
           ${SOURCES}
   )
   target_link_libraries(exodus_for_static PUBLIC Netcdf::all_libs)
   target_link_libraries(exodus_for_static PRIVATE SEACASExodus::all_libs)
   set_target_properties(exodus_for_static PROPERTIES OUTPUT_NAME exodus_for)
   # This keeps the library out of the `all_libs` targets...
   set_target_properties(exodus_for_static PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
   INSTALL(TARGETS exodus_for_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

   # For backwards compatibility, the Fortran exodus library was called libexoIIv2for.a
   # Create a symbolic link from libexodus_for.a to libexoIIv2for.a -- do only for static lib
   InstallSymLink(libexodus_for.a ${CMAKE_INSTALL_PREFIX}/lib/libexoIIv2for.a)
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
