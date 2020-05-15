INCLUDE(FortranSettings)

if (SEACASExodus_ENABLE_STATIC)
INCLUDE(InstallSymLink)
endif()

TRIBITS_SUBPACKAGE(Exodus_for)

INCLUDE_DIRECTORIES(
   "${CMAKE_CURRENT_SOURCE_DIR}/include"
   "${CMAKE_CURRENT_SOURCE_DIR}/test"
   )

SET(SOURCES "")
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

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
   set_target_properties(exodus_for_static PROPERTIES OUTPUT_NAME exodus_for)
   INSTALL(TARGETS exodus_for_static DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

   # For backwards compatibility, the Fortran exodus library was called libexoIIv2for.a
   # Create a symbolic link from libexodus_for.a to libexoIIv2for.a -- do only for static lib
   InstallSymLink(libexodus_for.a ${CMAKE_INSTALL_PREFIX}/lib/libexoIIv2for.a)
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
