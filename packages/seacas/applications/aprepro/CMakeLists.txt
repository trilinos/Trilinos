
TRIBITS_SUBPACKAGE(Aprepro)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

SET(SOURCES
	aprepro.C
)

TRIBITS_ADD_EXECUTABLE(
	aprepro
	NOEXEPREFIX
	NOEXESUFFIX
	SOURCES ${SOURCES}
)

install_executable(aprepro)

if(${CMAKE_PROJECT_NAME}_ENABLE_TESTS)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_unit_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro 
        ARGS -q ${CMAKE_CURRENT_SOURCE_DIR}/test.inp_app test.output
  OUTPUT_FILE test.stderr
  PASS_ANY
 TEST_1 CMND sed ARGS s+${CMAKE_CURRENT_SOURCE_DIR}++ 
		       ${CMAKE_CURRENT_BINARY_DIR}/test.stderr
  OUTPUT_FILE filtered.output
  PASS_ANY
 TEST_2 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test_standard.out
		       ${CMAKE_CURRENT_BINARY_DIR}/test.output
 TEST_3 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test.stderr.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/filtered.output
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_array_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro
	ARGS -q --include=${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test-array.i
	test-array.out
  OUTPUT_FILE test-array.stderr
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-array.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-array.out
 TEST_2 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-array.stderr.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-array.stderr
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_command_line_vars_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro
	ARGS -q var=42 svar='String' ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.i
	test-vars.out
  OUTPUT_FILE test-vars.stderr
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-vars.stderr
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_command_line_include_test
 TEST_0 CMND ${CMAKE_CURRENT_BINARY_DIR}/aprepro
	ARGS -q --include ${CMAKE_CURRENT_SOURCE_DIR}/test-include.i ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.i
	test-include.out
  OUTPUT_FILE test-include.stderr
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-include.stderr
 COMM serial		      
 FINAL_PASS_REGULAR_EXPRESSION 
)

endif()

TRIBITS_SUBPACKAGE_POSTPROCESS()
