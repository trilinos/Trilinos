INCLUDE(FortranSettings)

TRIBITS_SUBPACKAGE(ExoIIv2for32)

add_definitions(-DDEFAULT_REAL_INT)

INCLUDE_DIRECTORIES(
   "${CMAKE_CURRENT_SOURCE_DIR}/include"
   "${CMAKE_CURRENT_SOURCE_DIR}/test"
   )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

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
    set_target_properties(exoIIv2for32_static PROPERTIES OUTPUT_NAME exoIIv2for32)
    INSTALL(TARGETS exoIIv2for32_static DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
