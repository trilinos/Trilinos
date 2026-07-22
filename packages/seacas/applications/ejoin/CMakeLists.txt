
TRIBITS_SUBPACKAGE(Ejoin)

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.C)
TRIBITS_ADD_EXECUTABLE(
  ejoin
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES ${SOURCES}
  COMM serial mpi
  INSTALLABLE
  )

TRIBITS_SUBPACKAGE_POSTPROCESS()

if (${CMAKE_PROJECT_NAME}_ENABLE_SEACASExodiff)
TRIBITS_ADD_ADVANCED_TEST(ejoin_combine_nes_nodeset_merge
   TEST_0 EXEC ejoin ARGS --output combine_nes.g --offset 1,0,0 --match_nodeset_nodes p1:2,p2:1 --omit_nodeset p1:2,p2:1 --combine_element --combine_sideset --combine_nodeset ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC exodiff ARGS -show_all_diffs -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/2block-nes-gold.g combine_nes.g
       DIRECTORY ../exodiff
       NOEXEPREFIX NOEXESUFFIX
  COMM serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(ejoin_nodeset_merge
   TEST_0 EXEC ejoin ARGS --output ns-merge.g --offset p2:1,0,0 --match_nodeset_nodes p1:2,p2:1 --omit_nodeset p1:2,p2:1 ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC exodiff ARGS -show_all_diffs -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/2block-gold.g ns-merge.g
       DIRECTORY ../exodiff
       NOEXEPREFIX NOEXESUFFIX
  COMM serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(ejoin_combine_ne_nodeset_merge
   TEST_0 EXEC ejoin ARGS --output combine_ne.g --offset p1:0,0,0:p2:1,0,0 --match_nodeset_nodes p1:2,p2:1 --omit_nodeset p1:2,p2:1 --combine_element --combine_nodeset ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC exodiff ARGS -show_all_diffs -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/2block-ne-gold.g combine_ne.g
       DIRECTORY ../exodiff
       NOEXEPREFIX NOEXESUFFIX
  COMM serial
  XHOSTTYPE Windows
  )

TRIBITS_ADD_ADVANCED_TEST(ejoin_combine_ns_nodeset_merge
   TEST_0 EXEC ejoin ARGS -output combine_ns.g --offset 1,0,0 --match_nodeset_nodes p1:2,p2:1 --omit_nodeset p1:2,p2:1 --combine_sideset --combine_nodeset ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g ${CMAKE_CURRENT_SOURCE_DIR}/test/2block1.g
     NOEXEPREFIX NOEXESUFFIX
   TEST_1 EXEC exodiff ARGS -show_all_diffs -pedantic ${CMAKE_CURRENT_SOURCE_DIR}/test/2block-ns-gold.g combine_ns.g
       DIRECTORY ../exodiff
       NOEXEPREFIX NOEXESUFFIX
  COMM serial
  XHOSTTYPE Windows
  )
endif()
