INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(ExoIIv2for32)

add_definitions(-DDEFAULT_REAL_INT)

TRIBITS_INCLUDE_DIRECTORIES(
   "${CMAKE_CURRENT_SOURCE_DIR}/../exodus_for/include"
   "${CMAKE_CURRENT_SOURCE_DIR}/test"
   )

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_LIBRARY(
  exoIIv2for32
  HEADERS ""
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/exo_jack_32.c
  )

# Many users of the SEACAS installations still want a static exodus library...
if (SEACASExodus_ENABLE_STATIC)
    ADD_LIBRARY(
      exoIIv2for32_static
      STATIC
      ${CMAKE_CURRENT_SOURCE_DIR}/src/exo_jack_32.c
    )
    target_link_libraries(exoIIv2for32_static PUBLIC Netcdf::all_libs)
    target_link_libraries(exoIIv2for32_static PRIVATE SEACASExodus::all_libs)
    set_target_properties(exoIIv2for32_static PROPERTIES OUTPUT_NAME exoIIv2for32)
    # This keeps the library out of the `all_libs` targets...
    set_target_properties(exoIIv2for32_static PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
    INSTALL(TARGETS exoIIv2for32_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
