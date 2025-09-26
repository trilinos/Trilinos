# Copyright(C) 1999-2022, 2025 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details

TRIBITS_SUBPACKAGE(Aprepro)

TRIBITS_INCLUDE_DIRECTORIES(
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
	INSTALLABLE
)

# Currently doesn't run on windows due to use of sed...
TRIBITS_ADD_ADVANCED_TEST(
 aprepro_unit_test
 TEST_0 EXEC aprepro
        ARGS -q --legacy_output_format ${CMAKE_CURRENT_SOURCE_DIR}/test.inp_app test_unit.output
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test_standard.out
		       ${CMAKE_CURRENT_BINARY_DIR}/test_unit.output
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_array_test
 TEST_0 EXEC aprepro
	ARGS -q --legacy_output_format --info=aprepro-array.dump --include=${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test-array.i
	test-array.out
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-array.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-array.out
 TEST_2 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-array.stderr.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/aprepro-array.dump
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_command_line_vars_test
 TEST_0 EXEC aprepro
	ARGS -q var=42 svar='String' --info=test-vars.dump ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.i
	test-vars.out
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-vars.dump
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_command_line_include_test
 TEST_0 EXEC aprepro
	ARGS -q --info=test-include.dump --include ${CMAKE_CURRENT_SOURCE_DIR}/test-include.i ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.i
	test-include.out
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/test-dump.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test-include.dump
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

TRIBITS_ADD_ADVANCED_TEST(
 aprepro_test_dump_reread
 TEST_0 EXEC aprepro
        ARGS --info=test.dump -D ${CMAKE_CURRENT_SOURCE_DIR}/test.inp_app test_dump.output
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 EXEC aprepro
        ARGS --info=test-reread.dump -D test.dump test_dump.output
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_3 CMND diff ARGS -w test-filter.dump test-reread.dump
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

IF (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodus)
TRIBITS_ADD_ADVANCED_TEST(
 aprepro_test_exodus
 TEST_0 EXEC aprepro
        ARGS --legacy_output_format --include=${CMAKE_CURRENT_SOURCE_DIR} -q ${CMAKE_CURRENT_SOURCE_DIR}/exodus-test.apr test_exo.output
  NOEXEPREFIX NOEXESUFFIX
  PASS_ANY
 TEST_1 CMND diff ARGS -w
		       ${CMAKE_CURRENT_SOURCE_DIR}/exodus-test.gold
		       ${CMAKE_CURRENT_BINARY_DIR}/test_exo.output
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)
endif()

TRIBITS_SUBPACKAGE_POSTPROCESS()
