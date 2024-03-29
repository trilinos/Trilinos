
TRIBITS_SUBPACKAGE(Exodiff)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
SET(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ED_SystemInterface.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/exo_entity.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/exodiff.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/Tolerance.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/face_block.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/assembly.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/check.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/map.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/create_file.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/node_set.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/edge_block.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/side_set.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/exoII_read.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/stringx.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/exo_block.C
	    ${CMAKE_CURRENT_SOURCE_DIR}/util.C)

TRIBITS_ADD_EXECUTABLE(
  exodiff
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES ${SOURCES}
  COMM serial mpi
  INSTALLABLE
  )

TRIBITS_ADD_EXECUTABLE(
  test_generate_exo_files
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES test_generate_exo_files.c
  COMM serial mpi
  )

IF (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodus)

TRIBITS_ADD_ADVANCED_TEST(
 exodiff_test_edge_face
 TEST_0 EXEC ${CMAKE_CURRENT_BINARY_DIR}/test_generate_exo_files
        NOEXEPREFIX NOEXESUFFIX
 TEST_1 CMND ${CMAKE_CURRENT_BINARY_DIR}/exodiff
        ARGS ${CMAKE_CURRENT_BINARY_DIR}/edge-face-gold.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/edge-face-no-diff.exo
		    OUTPUT_FILE exodiff-no-diff.stdout
		    PASS_REGULAR_EXPRESSION "exodiff: Files are the same"
 TEST_2 CMND ${CMAKE_CURRENT_BINARY_DIR}/exodiff
        ARGS ${CMAKE_CURRENT_BINARY_DIR}/edge-face-gold.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/edge-diff.exo
		    OUTPUT_FILE exodiff-edge-diff.stdout
		    PASS_REGULAR_EXPRESSION "exodiff: Files are different"
 TEST_3 CMND ${CMAKE_CURRENT_BINARY_DIR}/exodiff
        ARGS ${CMAKE_CURRENT_BINARY_DIR}/edge-face-gold.exo
		         ${CMAKE_CURRENT_BINARY_DIR}/face-diff.exo
		    OUTPUT_FILE exodiff-face-diff.stdout
		    PASS_REGULAR_EXPRESSION "exodiff: Files are different"
 COMM mpi serial
 OVERALL_NUM_MPI_PROCS 1
 FINAL_PASS_REGULAR_EXPRESSION
 XHOSTTYPE Windows
)

endif()

TRIBITS_SUBPACKAGE_POSTPROCESS()
