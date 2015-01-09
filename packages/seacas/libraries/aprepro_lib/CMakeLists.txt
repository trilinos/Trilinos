
option(USE_TRIBITS "Use the Tribits modules to build Aprepro" ON)
if(USE_TRIBITS)
  INCLUDE(TribitsSubPackageMacros)
  INCLUDE(TribitsLibraryMacros)

  TRIBITS_SUBPACKAGE(Aprepro_lib)
endif()

SET(APREPRO_VERSION_MAJOR "4")
SET(APREPRO_VERSION_MINOR "02")
SET(APREPRO_VERSION_PATCH "0")
SET(APREPRO_VERSION "${APREPRO_VERSION_MAJOR}.${APREPRO_VERSION_MINOR}")
SET(APREPRO_VERSION_FULL "${APREPRO_VERSION}.${APREPRO_VERSION_PATCH}")

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

SET(HEADERS
  aprepro.h
)

SET(SOURCES
 apr_aprepro.cc apr_util.cc apr_stats.cc apr_init.cc apr_builtin.cc apr_units.cc apr_tokenize.cc
 apr_scanner.cc  apr_parser.cc apr_array.cc
)

option(GENERATE_FLEX_FILES "Generate source files using flex" OFF)
if(GENERATE_FLEX_FILES)
  find_program(FLEX_EXE NAMES flex
               DOC "Flex Executable")
  if(NOT FLEX_EXE)
    message(SEND_ERROR "The CMake variable FLEX_EXE does not point to a valid flex executable.")
  endif()

  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/apr_scanner.cc
                     COMMAND ${FLEX_EXE}
                     ARGS -o apr_scanner.cc
                          aprepro.ll
                     MAIN_DEPENDENCY aprepro.ll
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     )
endif()

option(GENERATE_BISON_FILES "Generate source files using bison" OFF)
if(GENERATE_BISON_FILES)
  find_program(BISON_EXE NAMES bison
               DOC "Bison Executable")
  if(NOT BISON_EXE)
    message(SEND_ERROR "The CMake variable BISON_EXE does not point to a valid bison executable.")
  endif()

  if(WIN32)
    get_filename_component(BISON_PATH ${BISON_EXE} PATH)
    set(ENV{PATH} "${BISON_PATH};$ENV{PATH}")
  endif(WIN32)

  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/aprepo_parser.h
                            ${CMAKE_CURRENT_SOURCE_DIR}/apr_parser.cc
                     COMMAND ${BISON_EXE}
                     ARGS -y --defines=aprepro_parser.h
                          -o apr_parser.cc
                          aprepro.yy
                     MAIN_DEPENDENCY aprepro.yy
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     )
endif()

if(USE_TRIBITS)
  TRIBITS_ADD_LIBRARY(
    aprepro_lib
    HEADERS	${HEADERS}
    SOURCES ${SOURCES}
  )
  TRIBITS_SUBPACKAGE_POSTPROCESS()

else()
  add_library(aprepro_lib ${HEADERS} ${SOURCES})
endif()

option(BUILD_TEST_APP "Build the test application" OFF)
if(BUILD_TEST_APP)
  add_executable(aprepro_test_app apr_test.cc)
  target_link_libraries(aprepro_test_app aprepro_lib)
endif()
