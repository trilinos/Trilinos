INCLUDE(TribitsAddExecutable)
SEACAS_PACKAGE(Aprepro)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
)

ADD_DEFINITIONS(-DNO_EXODUSII)

SET(HEADERS "")
SET(SOURCES
  ap_builtin.c
  ap_stat.c
  y.tab.c
  ap_util.c
  ap_init.c
  ap_hash.c
  ap_liballoc.c
  ap_units.c
  ap_aprepro.c
  ap_array.c
  )

TRIBITS_ADD_EXECUTABLE(
  aprepro
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES ${SOURCES}
  IMPORTEDLIBS m
  COMM serial mpi
  )
install_executable(aprepro)

option(GENERATE_FLEX_FILES "Generate source files using flex" OFF)
if(GENERATE_FLEX_FILES)
  find_program(FLEX_EXE NAMES flex
               DOC "Flex Executable")
  if(NOT FLEX_EXE)
    message(SEND_ERROR "The CMake variable FLEX_EXE does not point to a valid flex executable.")
  endif()

  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/lex.yy.c
                     COMMAND ${FLEX_EXE}
                     ARGS aprepro.l
                     MAIN_DEPENDENCY aprepro.l
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

  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/y.tab.h
                            ${CMAKE_CURRENT_SOURCE_DIR}/y.tab.c
                     COMMAND ${BISON_EXE}
                     ARGS --yacc --defines aprepro.y
                     MAIN_DEPENDENCY aprepro.y
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     )
endif()

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_unit_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro ARGS -q ${CMAKE_CURRENT_SOURCE_DIR}/test.inp_app
  OUTPUT_FILE test.output
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test_standard.out
		       ${CMAKE_CURRENT_BINARY_DIR}/test.output
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_array_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro 
	ARGS -q --include ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test-array.i
  OUTPUT_FILE test-array.out
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-array.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-array.out
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

SEACAS_POSTPROCESS()

